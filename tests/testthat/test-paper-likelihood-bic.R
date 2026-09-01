# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

test_that("observed likelihood is the log-sum-exp over latent states", {
  log_joint <- matrix(c(-2, -0.3, -1.1,
                        -0.4, -1.7, -0.8,
                        -3.2, -0.2, -2.1,
                        -0.7, -1.4, -0.5), nrow = 4, byrow = TRUE)
  expected <- sum(apply(log_joint, 1, function(x) log(sum(exp(x)))))
  r <- t(apply(log_joint, 1, lse_vec))

  expect_equal(observed_loglik(log_joint), expected, tolerance = 1e-12)
  expect_equal(expected - sum(log_joint * r), -sum(r * log(r)), tolerance = 1e-12)
})

test_that("parallel observed likelihood enumerates every joint state", {
  log_joint <- array(seq(-3, 1, length.out = 2 * 3 * 5), dim = c(2, 3, 5))
  expected <- sum(vapply(1:5, function(i) {
    x <- as.numeric(log_joint[, , i])
    max(x) + log(sum(exp(x - max(x))))
  }, numeric(1)))

  expect_equal(
    observed_loglik(log_joint, observations = "last"),
    expected,
    tolerance = 1e-12
  )
  expect_equal(cal_loglik(log_joint), expected, tolerance = 1e-12)
})

test_that("one-layer parallel likelihood uses columns as observations", {
  log_joint <- matrix(c(log(0.1), log(0.3), log(0.6),
                        log(0.7), log(0.2), log(0.1)), nrow = 3)
  expected <- sum(apply(log_joint, 2, LUCIDus:::safe_log_sum_exp))

  expect_equal(cal_loglik(log_joint), expected, tolerance = 1e-12)
})

test_that("early fit reports likelihood from its final relabeled parameters", {
  set.seed(9301)
  n <- 180
  G <- cbind(g1 = rnorm(n), g2 = rnorm(n))
  X <- rbinom(n, 1, plogis(-0.2 + G[, 1])) + 1L
  Z <- cbind(z1 = rnorm(n, c(-1.5, 1.5)[X], 0.6),
             z2 = rnorm(n, c(-1, 1)[X], 0.7))
  age <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "age"))
  Y <- c(-0.6, 1.1)[X] + 0.4 * age[, 1] + rnorm(n, sd = 0.7)
  fit <- NULL
  invisible(capture.output(fit <- est_lucid(
    "early", G, Z, Y, CoY = age, K = 2, family = "normal",
    init_omic.data.model = "VVV", seed = 21, max_itr = 80,
    max_tot.itr = 160
  )))
  sigma <- array(unlist(fit$res_Sigma), dim = c(ncol(Z), ncol(Z), 2))
  final_joint <- Estep_early(
    beta = fit$res_Beta, mu = fit$res_Mu, sigma = sigma,
    gamma = fit$res_Gamma, G = G, Z = fit$Z, Y = matrix(Y), CoY = age,
    N = n, K = 2, family.list = normal(2, 1), itr = 2,
    useY = TRUE, dimCoY = 1, ind.na = check_na(fit$Z)$indicator_na
  )

  expect_equal(fit$likelihood, observed_loglik(final_joint), tolerance = 1e-8)
  expect_equal(fit$inclusion.p, t(apply(final_joint, 1, lse_vec)),
               tolerance = 1e-8)
})

test_that("parallel fit reports likelihood from its final relabeled parameters", {
  set.seed(9302)
  n <- 200
  G <- cbind(g1 = rnorm(n), g2 = rnorm(n))
  X1 <- rbinom(n, 1, plogis(-0.3 + G[, 1])) + 1L
  X2 <- rbinom(n, 1, plogis(0.2 - G[, 2])) + 1L
  Z <- list(
    cbind(z11 = rnorm(n, c(-1.4, 1.4)[X1], 0.6),
          z12 = rnorm(n, c(-0.8, 0.8)[X1], 0.7)),
    cbind(z21 = rnorm(n, c(-1.3, 1.3)[X2], 0.6),
          z22 = rnorm(n, c(-0.9, 0.9)[X2], 0.7))
  )
  age <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "age"))
  Y <- -0.5 + 0.8 * (X1 == 2) + 0.6 * (X2 == 2) +
    0.3 * age[, 1] + rnorm(n, sd = 0.75)
  fit <- NULL
  invisible(capture.output(fit <- est_lucid(
    "parallel", G, Z, Y, CoY = age, K = c(2, 2), family = "normal",
    init_omic.data.model = "VVV", seed = 22, max_itr = 80,
    max_tot.itr = 160
  )))
  final_joint <- Estep(
    G = G, Z = fit$Z, Y = matrix(Y), Beta = fit$res_Beta$Beta,
    Mu = fit$res_Mu, Sigma = fit$res_Sigma,
    Delta = fit$res_Gamma$Gamma, family = "normal", useY = TRUE,
    na_pattern = lapply(fit$Z, check_na), CoY = age
  )
  final_r <- Estep_to_r(final_joint, K = c(2, 2), N = n)

  expect_equal(fit$likelihood, cal_loglik(final_joint), tolerance = 1e-8)
  expect_equal(fit$z, final_r, tolerance = 1e-8)
})

test_that("early BIC counts symmetric covariance and adjustment covariates", {
  K <- 3
  object <- structure(list(
    select = list(selectG = c(TRUE, FALSE), selectZ = c(TRUE, TRUE, FALSE)),
    K = K,
    res_Gamma = early_gamma_from_levels(c(-1, 0, 1), covariate = 0.2,
                                         sigma = c(1, 1.2, 0.8),
                                         covariate_names = "age"),
    res_Beta = matrix(0, K, 4),
    res_Mu = matrix(0, K, 3),
    res_Sigma = replicate(K, diag(3), simplify = FALSE),
    likelihood = -50,
    inclusion.p = matrix(1 / K, 20, K),
    family = "normal", useY = TRUE,
    Z = matrix(0, 20, 3),
    var.names = list(Gnames = c("g1", "g2", "sex"), Znames = paste0("z", 1:3)),
    missing_summary = NULL, Rho = NULL
  ), class = "early_lucid")
  # Eq 18: BICp = -2logL + (D - DG - DZ) log N, where D is the FULL model
  # dimension and DG / DZ count the variables whose effects are zero in every
  # cluster.  Eq 18 subtracts one per deselected VARIABLE; it does not delete
  # that variable's entire mean/covariance/regression block, which is what the
  # previous implementation did (and what this expectation used to encode).
  P <- 2L; nCoG <- 1L; M <- 3L
  D <- (P + nCoG + 1L) * (K - 1L) + K * M + K * M * (M + 1L) / 2 + (K + 1L + K)
  DG <- 1L  # selectG = c(TRUE, FALSE)
  DZ <- 1L  # selectZ = c(TRUE, TRUE, FALSE)
  expected <- D - DG - DZ

  expect_equal(early_npars(object), expected)
  expect_equal(summary_lucid(object)$model_fit$n_parameters, expected)
})

test_that("parallel BIC handles unequal dimensions and layer-specific selection", {
  K <- c(2, 3)
  Delta <- parallel_delta_from_coef(c(-0.5, 0.4, -0.2, 0.7, 0.3), K,
                                    "binomial", covariate_names = "age")
  object <- structure(list(
    K = K, N = 30,
    select = list(
      selectG = c(TRUE, TRUE),
      selectG_layer = list(c(TRUE, FALSE), c(TRUE, TRUE)),
      selectZ = list(matrix(TRUE, 2, 2), matrix(TRUE, 3, 1))
    ),
    var.names = list(Gnames = c("g1", "g2", "sex"),
                     Znames = list(c("z11", "z12"), "z21")),
    Z = list(matrix(0, 30, 2), matrix(0, 30, 1)),
    res_Gamma = list(Gamma = Delta, fit = Delta$fit),
    inclusion.p = list(matrix(0.5, 30, 2), matrix(1 / 3, 30, 3)),
    likelihood = -80, family = "binomial", useY = TRUE
  ), class = "lucid_parallel")
  expected_g <- (1 + 1 + 1) * 1 + (2 + 1 + 1) * 2
  expected_z <- 2 * 2 + 2 * 2 * 3 / 2 + 3 * 1 + 3 * 1 * 2 / 2
  expected <- expected_g + expected_z + length(Delta$coef)

  expect_equal(parallel_npars(object), expected)
  expect_equal(cal_bic_parallel(object), -2 * object$likelihood + expected * log(30))
})

test_that("unsupervised parameter counts exclude post-hoc outcome coefficients", {
  object <- structure(list(
    select = list(selectG = TRUE, selectZ = TRUE), K = 2,
    res_Gamma = early_gamma_from_levels(c(-1, 1)),
    res_Beta = matrix(0, 2, 2), res_Mu = matrix(0, 2, 1),
    res_Sigma = list(matrix(1), matrix(1)), likelihood = -10,
    inclusion.p = matrix(.5, 10, 2), family = "binary", useY = FALSE,
    Z = matrix(0, 10, 1), var.names = list(Gnames = "G", Znames = "Z"),
    missing_summary = NULL, Rho = NULL
  ), class = "early_lucid")

  expect_equal(early_npars(object), 2 + 2 + 2)
})

test_that("unsupervised early fitting is independent of Y before final risk ordering", {
  set.seed(201)
  n <- 70
  G <- matrix(rnorm(n * 2), nrow = n)
  Z <- matrix(rnorm(n * 4), nrow = n)
  Y <- rbinom(n, 1, 0.4)
  args <- list(lucid_model = "early", G = G, Z = Z, K = 2, useY = FALSE,
               family = "binary", init_par = "random", seed = 44,
               max_itr = 8, max_tot.itr = 16, relabel_by_outcome = FALSE)
  invisible(capture.output(a <- do.call(est_lucid, c(args, list(Y = Y)))))
  invisible(capture.output(b <- do.call(est_lucid, c(args, list(Y = rev(Y))))))

  expect_equal(a$res_Beta, b$res_Beta)
  expect_equal(a$res_Mu, b$res_Mu)
  expect_equal(a$res_Sigma, b$res_Sigma)
  expect_equal(a$inclusion.p, b$inclusion.p)
})

test_that("early binary summary reports baseline odds and the LC odds ratio", {
  gamma <- early_gamma_from_levels(c(-1.5, 0.2))
  text <- capture.output(f.binary.early(gamma, 2, NULL))
  expect_true(any(grepl(sprintf("%.4f", exp(-1.5)), text, fixed = TRUE)))
  expect_equal(unname(gamma$beta[2]), 1.7)
  expect_equal(unname(exp(gamma$beta[2])), exp(0.2) / exp(-1.5))
})

test_that("unsupervised early binary report uses weighted latent-state regression", {
  set.seed(202)
  n <- 120
  latent <- rep(0:1, each = n / 2)
  G <- cbind(rnorm(n, latent), rnorm(n))
  Z <- cbind(rnorm(n, 2 * latent, .6), rnorm(n, -latent, .6), rnorm(n))
  CoY <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "age"))
  Y <- rbinom(n, 1, plogis(-1 + 1.1 * latent + .35 * CoY[, 1]))
  invisible(capture.output(fit <- estimate_lucid(
    lucid_model = "early", G = G, Z = Z, Y = Y, CoY = CoY, K = 2,
    useY = FALSE, family = "binary", init_par = "random", seed = 73,
    max_itr = 20, max_tot.itr = 40
  )))
  long <- data.frame(Y = rep(Y, each = 2),
                     state = factor(rep(1:2, times = n), levels = 1:2),
                     age = rep(CoY[, 1], each = 2),
                     w = as.vector(t(fit$inclusion.p)))
  expected <- suppressWarnings(glm(Y ~ state + age, family = binomial(),
                                    data = long, weights = w))

  expect_equal(unname(fit$res_Gamma$beta), unname(coef(expected)), tolerance = 1e-8)
  expect_equal(unname(exp(fit$res_Gamma$beta[2])),
               unname(exp(coef(expected)[2])), tolerance = 1e-8)
})

test_that("parallel outcome prediction integrates state-specific responses", {
  K <- c(2, 2)
  Delta <- parallel_delta_from_coef(c(-1, .8, -.4), K, "binomial")
  r <- array(c(.1, .2, .3, .4, .4, .3, .2, .1), dim = c(K, 2))
  expected <- colSums(matrix(r, 4, 2) * matrix(plogis(parallel_state_eta(Delta, N = 2)), 4, 2))
  weighted_logit <- plogis(colSums(matrix(r, 4, 2) * matrix(parallel_state_eta(Delta, N = 2), 4, 2)))

  expect_false(isTRUE(all.equal(expected, weighted_logit)))
  expect_true(all(expected > 0 & expected < 1))
})

test_that("unsupervised parallel report uses weighted joint-state regression", {
  set.seed(203)
  n <- 90
  G <- matrix(rnorm(n * 2), nrow = n)
  Z <- list(matrix(rnorm(n * 4), nrow = n), matrix(rnorm(n * 4), nrow = n))
  CoY <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "age"))
  Y <- rbinom(n, 1, plogis(-.6 + .3 * G[, 1] + .25 * CoY[, 1]))
  invisible(capture.output(fit <- estimate_lucid(
    lucid_model = "parallel", G = G, Z = Z, Y = Y, CoY = CoY,
    K = c(2, 2), useY = FALSE, family = "binary", init_par = "random",
    seed = 74, max_itr = 12, max_tot.itr = 24
  )))
  grid <- expand.grid(1:2, 1:2)
  long <- data.frame(
    Y = rep(Y, each = 4),
    Layer1_LC2 = rep(as.numeric(grid[[1]] == 2), times = n),
    Layer2_LC2 = rep(as.numeric(grid[[2]] == 2), times = n),
    age = rep(CoY[, 1], each = 4),
    w = as.vector(matrix(fit$z, 4, n))
  )
  expected <- suppressWarnings(glm(Y ~ Layer1_LC2 + Layer2_LC2 + age,
                                    data = long, weights = w, family = binomial()))

  expect_equal(unname(parallel_delta_coef(fit$res_Gamma$Gamma)),
               unname(coef(expected)), tolerance = 1e-8)
})

test_that("unsupervised BIC excludes post-hoc outcome parameters", {
  object <- structure(list(
    select = list(selectG = TRUE, selectZ = TRUE), K = 2,
    res_Gamma = early_gamma_from_levels(c(-1, 1)),
    res_Beta = matrix(0, 2, 2), res_Mu = matrix(0, 2, 1),
    res_Sigma = list(matrix(1), matrix(1)), likelihood = -10,
    inclusion.p = matrix(.5, 10, 2), family = "binary", useY = FALSE,
    Z = matrix(0, 10, 1), var.names = list(Gnames = "G", Znames = "Z"),
    missing_summary = NULL, Rho = NULL
  ), class = "early_lucid")
  out <- summary_lucid(object)
  expected_without_y <- (1 + 1) * (2 - 1) + (1 * 2 + 1^2 * 2)
  expect_equal(out$model_fit$n_parameters, expected_without_y)
})

test_that("serial upstream stages are not outcome-relabelled", {
  set.seed(204)
  n <- 55
  G <- matrix(rnorm(n * 2), nrow = n)
  Z <- list(matrix(rnorm(n * 3), nrow = n), matrix(rnorm(n * 3), nrow = n))
  Y <- rnorm(n)
  args <- list(lucid_model = "serial", G = G, Z = Z, K = list(2, 2),
               useY = FALSE, family = "normal", init_par = "random",
               seed = 81, max_itr = 6, max_tot.itr = 12)
  invisible(capture.output(a <- do.call(estimate_lucid, c(args, list(Y = Y)))))
  invisible(capture.output(b <- do.call(estimate_lucid, c(args, list(Y = rev(Y))))))

  expect_equal(a$submodel[[1]]$res_Beta, b$submodel[[1]]$res_Beta)
  expect_equal(a$submodel[[1]]$res_Mu, b$submodel[[1]]$res_Mu)
  expect_equal(a$submodel[[1]]$inclusion.p, b$submodel[[1]]$inclusion.p)
  expect_false(identical(a$submodel[[1]]$res_Gamma$beta,
                         b$submodel[[1]]$res_Gamma$beta))
})

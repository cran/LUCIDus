# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# Regression tests for critical parallel-model paths

test_that("summary_lucid works with current parallel select shapes", {
  set.seed(1008)
  G <- matrix(rnorm(240), nrow = 60)
  Z1 <- matrix(rnorm(600), nrow = 60)
  Z2 <- matrix(rnorm(600), nrow = 60)
  Z <- list(Z1 = Z1, Z2 = Z2)
  Y <- rnorm(60)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      seed = 1008,
      useY = TRUE
    )
  )))

  s <- summary_lucid(fit)
  expect_s3_class(fit, "lucid_parallel")
  expect_s3_class(s, "sumlucid_parallel")
  expect_true(is.finite(s$model_fit$BIC))
})

test_that("summary_lucid parallel fallback works when selectG is NULL", {
  set.seed(1008)
  G <- matrix(rnorm(240), nrow = 60)
  Z1 <- matrix(rnorm(600), nrow = 60)
  Z2 <- matrix(rnorm(600), nrow = 60)
  Z <- list(Z1 = Z1, Z2 = Z2)
  Y <- rnorm(60)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      seed = 1008,
      useY = TRUE
    )
  )))

  fit$select$selectG <- NULL
  fit$select$selectG_layer <- NULL
  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_parallel")
  expect_equal(s$model_info$n_features$G, length(fit$var.names$Gnames))
})

test_that("parallel prediction for 2 layers matches manual state integration", {
  set.seed(1008)
  G <- matrix(rnorm(240), nrow = 60)
  Z1 <- matrix(rnorm(600), nrow = 60)
  Z2 <- matrix(rnorm(600), nrow = 60)
  Z <- list(Z1 = Z1, Z2 = Z2)
  Y <- rnorm(60)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      seed = 1008,
      useY = TRUE
    )
  )))

  pred <- predict_lucid(
    model = fit,
    lucid_model = "parallel",
    G = G,
    Z = Z,
    Y = Y,
    response = FALSE
  )

  na_pattern <- lapply(Z, check_na)
  likelihood <- Estep(G, Z, matrix(Y, ncol = 1), fit$res_Beta$Beta,
                      fit$res_Mu, fit$res_Sigma, fit$res_Gamma$Gamma,
                      fit$family, TRUE, na_pattern)
  r <- Estep_to_r(likelihood, fit$K, nrow(G))
  manual_y <- colSums(matrix(r, nrow = prod(fit$K), ncol = nrow(G)) *
                      matrix(parallel_state_eta(fit$res_Gamma$Gamma, N = nrow(G)),
                             nrow = prod(fit$K), ncol = nrow(G)))

  expect_equal(as.vector(pred$pred.y), manual_y, tolerance = 1e-7)
})

test_that("parallel E-step remains finite with all-missing rows in one layer", {
  skip_if_not_installed("mix")
  set.seed(1008)
  G <- matrix(rnorm(320), nrow = 80)
  Z1 <- matrix(rnorm(800), nrow = 80)
  Z2 <- matrix(rnorm(800), nrow = 80)
  Z1[1:3, ] <- NA
  Z1[4:8, 1:2] <- NA
  Z <- list(Z1 = Z1, Z2 = Z2)
  Y <- rnorm(80)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      seed = 1008,
      init_impute = "mix",
      useY = TRUE
    )
  )))

  expect_s3_class(fit, "lucid_parallel")
  expect_equal(length(fit$inclusion.p), 2)
  for (i in 1:2) {
    expect_true(all(is.finite(fit$inclusion.p[[i]])))
    expect_equal(rowSums(fit$inclusion.p[[i]]), rep(1, nrow(G)), tolerance = 1e-6)
  }
})

test_that("parallel missing-data path keeps all-missing rows as NA and imputes partial rows", {
  skip_if_not_installed("mix")
  set.seed(1008)
  G <- matrix(rnorm(240), nrow = 60)
  Z1 <- matrix(rnorm(600), nrow = 60)
  Z2 <- matrix(rnorm(600), nrow = 60)
  Z1[1, ] <- NA        # pattern 3 (all missing)
  Z1[2, 1:3] <- NA     # pattern 2 (partial missing)
  Z <- list(Z1 = Z1, Z2 = Z2)
  Y <- rnorm(60)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      seed = 1008,
      init_impute = "mix",
      useY = TRUE
    )
  )))

  expect_true(all(is.na(fit$Z[[1]][1, ])))
  expect_true(all(is.finite(fit$Z[[1]][2, ])))
})

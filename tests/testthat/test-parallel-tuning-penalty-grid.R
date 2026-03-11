# Parallel tuning tests for K and regularization grids

test_that("tune_lucid parallel supports regularization grids", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(240), nrow = 40)
  Z2 <- matrix(rnorm(240), nrow = 40)
  Z <- list(Z1, Z2)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    tuned <- tune_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = list(2:3, 2),
      Rho_G = c(0, 0.05),
      Rho_Z_Mu = c(0, 0.1),
      Rho_Z_Cov = 0,
      max_itr = 8,
      tol = 1e-1,
      seed = 1008,
      useY = TRUE
    )
  )))

  expect_true(is.data.frame(tuned$tune_K))
  expect_equal(nrow(tuned$tune_K), 8)
  expect_true(all(c("K1", "K2", "Rho_G", "Rho_Z_Mu", "Rho_Z_Cov", "BIC") %in% colnames(tuned$tune_K)))
  expect_s3_class(tuned$model_opt, "lucid_parallel")
})

test_that("lucid wrapper parallel can tune over penalty vectors", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(240), nrow = 40)
  Z2 <- matrix(rnorm(240), nrow = 40)
  Z <- list(Z1, Z2)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = list(2:3, 2),
      Rho_G = c(0, 0.05),
      Rho_Z_Mu = 0,
      Rho_Z_Cov = 0,
      max_itr = 8,
      tol = 1e-1,
      seed = 1008,
      useY = TRUE
    )
  )))

  expect_s3_class(fit, "lucid_parallel")
  expect_true(fit$K[1] %in% c(2, 3))
  expect_equal(fit$K[2], 2)
  expect_true(fit$Rho$Rho_G %in% c(0, 0.05))
})

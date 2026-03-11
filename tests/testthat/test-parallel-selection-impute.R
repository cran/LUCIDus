test_that("parallel LUCID feature selection works", {
  set.seed(123)
  N <- 30

  G <- matrix(rnorm(N * 2), nrow = N)
  Z1 <- matrix(rnorm(N * 3), nrow = N)
  Z2 <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      family = "normal",
      K = c(2, 2),
      init_impute = "mix",
      Rho_G = 0.1,
      Rho_Z_Mu = 0.1,
      Rho_Z_Cov = 0,
      max_itr = 6,
      tol = 1e-1,
      seed = 123
    )
  )))

  expect_s3_class(fit, "lucid_parallel")
  expect_type(fit$select$selectG, "logical")
  expect_equal(length(fit$select$selectG), ncol(G))
  expect_true(is.list(fit$select$selectG_layer))
  expect_true(is.list(fit$select$selectZ))
  expect_equal(length(fit$select$selectG_layer), 2)
  expect_equal(length(fit$select$selectZ), 2)

  expect_true(all(sapply(fit$res_Beta$Beta, function(x) all(is.finite(x)))))
  expect_true(all(sapply(fit$res_Mu, function(x) all(is.finite(x)))))
  expect_true(all(sapply(fit$res_Sigma, function(x) all(is.finite(x)))))
})

test_that("parallel LUCID missing value imputation works", {
  set.seed(123)
  N <- 30

  G <- matrix(rnorm(N * 2), nrow = N)
  Z1 <- matrix(rnorm(N * 3), nrow = N)
  Z2 <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  # mixed missingness
  Z1[1, ] <- NA
  Z1[2, 1] <- NA
  Z2[3, ] <- NA
  Z2[4, 2] <- NA

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      family = "normal",
      K = c(2, 2),
      init_impute = "mix",
      max_itr = 6,
      tol = 1e-1,
      seed = 123
    )
  )))

  # Under mix policy, all-missing rows remain NA while partial-missing rows are imputed.
  expect_true(all(is.na(fit$Z[[1]][1, ])))
  expect_true(all(is.na(fit$Z[[2]][3, ])))
  expect_true(is.finite(fit$Z[[1]][2, 1]))
  expect_true(is.finite(fit$Z[[2]][4, 2]))
})

test_that("parallel Z-penalty selection works with listwise-missing rows", {
  set.seed(321)
  N <- 40

  G <- matrix(rnorm(N * 2), nrow = N)
  Z1 <- matrix(rnorm(N * 4), nrow = N)
  Z2 <- matrix(rnorm(N * 4), nrow = N)
  Y <- rnorm(N)

  # Include both listwise and sporadic missingness.
  Z1[1, ] <- NA
  Z1[2, 1] <- NA
  Z2[3, ] <- NA
  Z2[4, 2] <- NA

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      family = "normal",
      K = c(2, 2),
      init_impute = "mix",
      Rho_G = 0,
      Rho_Z_Mu = 0.1,
      Rho_Z_Cov = 0.05,
      max_itr = 8,
      tol = 1e-1,
      seed = 321
    )
  )))

  expect_s3_class(fit, "lucid_parallel")
  expect_true(all(sapply(fit$res_Mu, function(x) all(is.finite(x)))))
  expect_true(all(sapply(fit$res_Sigma, function(x) all(is.finite(x)))))
  expect_true(is.list(fit$select$selectZ))
})

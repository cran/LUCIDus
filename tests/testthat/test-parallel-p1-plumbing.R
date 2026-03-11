# Regression tests for parallel-model parameter plumbing

test_that("parallel model records and runs requested init_par", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z <- list(matrix(rnorm(400), nrow = 40), matrix(rnorm(400), nrow = 40))
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit_random <- est_lucid(
      lucid_model = "parallel",
      G = G, Z = Z, Y = Y,
      family = "normal",
      K = c(2, 2),
      init_par = "random",
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    fit_mclust <- est_lucid(
      lucid_model = "parallel",
      G = G, Z = Z, Y = Y,
      family = "normal",
      K = c(2, 2),
      init_par = "mclust",
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  expect_equal(fit_random$init_par, "random")
  expect_equal(fit_mclust$init_par, "mclust")
})

test_that("parallel model stores EM control settings for bootstrap reuse", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z <- list(matrix(rnorm(400), nrow = 40), matrix(rnorm(400), nrow = 40))
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- est_lucid(
      lucid_model = "parallel",
      G = G, Z = Z, Y = Y,
      family = "normal",
      K = c(2, 2),
      tol = 1e-2,
      max_itr = 7,
      max_tot.itr = 25,
      seed = 1008
    )
  )))

  expect_true(is.list(fit$em_control))
  expect_equal(fit$em_control$tol, 1e-2)
  expect_equal(fit$em_control$max_itr, 7)
  expect_equal(fit$em_control$max_tot.itr, 25)
})

test_that("parallel Z penalties are plumbed into fit metadata and select structure", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(400), nrow = 40)
  Z2 <- matrix(rnorm(400), nrow = 40)
  Z <- list(Z1, Z2)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- est_lucid(
      lucid_model = "parallel",
      G = G, Z = Z, Y = Y,
      family = "normal",
      K = c(2, 2),
      Rho_Z_Mu = 0.2,
      Rho_Z_Cov = 0.05,
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  expect_equal(fit$Rho$Rho_Z_Mu, 0.2)
  expect_equal(fit$Rho$Rho_Z_Cov, 0.05)
  expect_true(is.list(fit$select$selectZ))
  expect_equal(length(fit$select$selectZ), 2)
  expect_true(all(sapply(fit$select$selectZ, is.matrix)))
  expect_equal(dim(fit$select$selectZ[[1]]), c(2, ncol(Z1)))
  expect_equal(dim(fit$select$selectZ[[2]]), c(2, ncol(Z2)))
})

test_that("parallel G penalties return overall and per-layer selectG objects", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z <- list(matrix(rnorm(400), nrow = 40), matrix(rnorm(400), nrow = 40))
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- est_lucid(
      lucid_model = "parallel",
      G = G, Z = Z, Y = Y,
      family = "normal",
      K = c(2, 2),
      Rho_G = 0.1,
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  expect_type(fit$select$selectG, "logical")
  expect_equal(length(fit$select$selectG), ncol(G))
  expect_true(is.list(fit$select$selectG_layer))
  expect_equal(length(fit$select$selectG_layer), 2)
  expect_true(all(sapply(fit$select$selectG_layer, is.logical)))
  expect_equal(length(fit$select$selectG_layer[[1]]), ncol(G))
  expect_equal(length(fit$select$selectG_layer[[2]]), ncol(G))
  expect_equal(fit$select$selectG,
               fit$select$selectG_layer[[1]] | fit$select$selectG_layer[[2]])
})

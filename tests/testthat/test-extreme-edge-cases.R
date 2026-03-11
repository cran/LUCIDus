test_that("early model handles p >> n exposure setting with selection", {
  set.seed(77)
  N <- 20
  G <- matrix(rnorm(N * 60), nrow = N)  # p >> n
  Z <- matrix(rnorm(N * 6), nrow = N)
  Y <- rnorm(N)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "early",
      G = G,
      Z = Z,
      Y = Y,
      family = "normal",
      K = 2,
      Rho_G = 0.05,
      init_omic.data.model = NULL,
      max_itr = 4,
      tol = 1e-1,
      seed = 77
    )
  )))

  expect_s3_class(fit, "early_lucid")
  expect_true(is.finite(fit$likelihood))
  expect_equal(length(fit$select$selectG), ncol(G))
  expect_true(any(fit$select$selectG))
})

test_that("early binary model remains stable under near-separation", {
  set.seed(88)
  N <- 80
  G <- matrix(rnorm(N * 6), nrow = N)
  Z <- matrix(rnorm(N * 8), nrow = N)
  lin <- 6 * scale(G[, 1]) + 0.2 * rnorm(N)
  Y <- as.numeric(lin > 0)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "early",
      G = G,
      Z = Z,
      Y = Y,
      family = "binary",
      K = 2,
      init_omic.data.model = NULL,
      max_itr = 8,
      tol = 1e-1,
      seed = 88
    )
  )))

  expect_s3_class(fit, "early_lucid")
  expect_true(is.finite(fit$likelihood))
  expect_true(all(is.finite(fit$res_Gamma$beta)))
  expect_true(all(is.finite(fit$inclusion.p)))
  expect_equal(nrow(fit$inclusion.p), N)
})

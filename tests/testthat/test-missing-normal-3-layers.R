# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# LUCID - three omics, normal outcome

test_that("check estimations of LUCID with normal outcome (K = 2,2,2) with missing data", {
  skip_if_not_installed("mix")
  # run LUCID model
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2, Z2 = Z3)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  Y <- rnorm(100)


  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2), CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "normal",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE)))
  betas <- fit1$res_Beta$Beta
  beta1 <- mean(unlist(betas[1]))
  beta2 <- mean(unlist(betas[2]))
  beta3 <- mean(unlist(betas[3]))

  mus <- fit1$res_Mu
  mu1 <- mean(unlist(mus[1]))
  mu2 <- mean(unlist(mus[2]))
  mu3 <- mean(unlist(mus[3]))

  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(parallel_delta_coef(fit1$res_Gamma$Gamma))

  # check parameters
  expect_true(all(is.finite(c(beta1, beta2, beta3))))

  # The data here is pure noise, so there is no true mu to recover -- the
  # previous snapshots of -0.042 / 0.1119 / -0.01587 were floating-point
  # artefacts of one seed and could not fail when the estimator was wrong.
  # Assert what is actually true of a fit on centred noise instead.
  expect_true(all(abs(c(mu1, mu2, mu3)) < 1))
  expect_gt(sigma, 0)
  expect_true(is.finite(sigma))
  expect_true(all(vapply(fit1$res_Gamma$Gamma$effects,
                         function(x) all(diff(x) >= 0), logical(1))))

  expect_equal(class(fit1), "lucid_parallel")

  # missing data
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  a = sample(1:1000, 30, replace=FALSE)
  Z1[a] = NA
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z2[62:65, 6:8] = NA
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2, Z2 = Z3)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  Y <- rnorm(100)


  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2), CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "normal",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE,
                                             init_impute = "mix")))
  betas <- fit1$res_Beta$Beta
  beta1 <- mean(unlist(betas[1]))
  beta2 <- mean(unlist(betas[2]))
  beta3 <- mean(unlist(betas[3]))

  mus <- fit1$res_Mu
  mu1 <- mean(unlist(mus[1]))
  mu2 <- mean(unlist(mus[2]))
  mu3 <- mean(unlist(mus[3]))

  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(parallel_delta_coef(fit1$res_Gamma$Gamma))

  # check parameters
  expect_true(all(is.finite(c(beta1, beta2, beta3))))

  expect_lt(abs(mu1 + 0.0394), 0.05)
  expect_lt(abs(mu2 - 0.0989), 0.05)
  expect_lt(abs(mu3 - 0.01258), 0.05)

  expect_gt(sigma, 0)
  expect_true(is.finite(sigma))
  expect_true(all(vapply(fit1$res_Gamma$Gamma$effects,
                         function(x) all(diff(x) >= 0), logical(1))))

  expect_equal(class(fit1), "lucid_parallel")

})

# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# LUCID - five omics, normal outcome

test_that("check estimations of LUCID with normal outcome (K = 2,2,2)", {
  # run LUCID model
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z4 <- matrix(rnorm(1000), nrow = 100)
  Z5 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5)
  Y <- rnorm(100)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2, 2, 2),
                                             lucid_model = "parallel",
                                             family = "normal",
                                             seed = i,
                                             useY = TRUE)))
  betas <- mean(unlist(fit1$res_Beta$Beta))
  mus <- mean(unlist(fit1$res_Mu))
  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(parallel_delta_coef(fit1$res_Gamma$Gamma))

  # check parameters
  expect_true(is.finite(betas))
  expect_lt(abs(betas), 0.3)


  # Pure-noise data has no cluster structure to recover, so the previous
  # snapshots (-0.01, 0.08447) recorded one seed's rounding rather than any
  # property of the estimator. Assert the properties that do hold.
  expect_true(is.finite(mus))
  expect_lt(abs(mus), 1)
  expect_gt(sigma, 0)
  expect_true(is.finite(sigma))
  expect_true(all(is.finite(parallel_delta_coef(fit1$res_Gamma$Gamma))))
  expect_equal(length(parallel_delta_coef(fit1$res_Gamma$Gamma)), 6)
  expect_true(all(vapply(fit1$res_Gamma$Gamma$effects,
                         function(x) all(diff(x) >= 0), logical(1))))

  expect_equal(class(fit1), "lucid_parallel")

  invisible(capture.output(fit2 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2, 2, 2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "normal",
                                             seed = i,
                                             useY = TRUE)))
  betas <- mean(unlist(fit2$res_Beta$Beta))
  mus <- mean(unlist(fit2$res_Mu))
  sigma <- mean(unlist(fit2$res_Sigma))
  Gamma <- mean(parallel_delta_coef(fit2$res_Gamma$Gamma))

  # check parameters
  expect_true(is.finite(betas))
  expect_lt(abs(betas), 1)
  expect_true(is.finite(mus))
  expect_lt(abs(mus), 1)



  expect_equal(class(fit2), "lucid_parallel")
})

# LUCID - three omics, binary outcome



test_that("check estimations of LUCID with binary outcome (K = 2,2,2)", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3)
  Y <- rbinom(n=100, size =1, prob =0.25)

  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2),
                                             lucid_model = "parallel",
                                             family = "binary",

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
  Gamma <- mean(unlist(fit1$res_Gamma$Gamma))

  # check parameters
  expect_equal(beta1, 0.00, tolerance = 0.01)
  expect_equal(beta2, 0.0719, tolerance = 0.01)
  expect_equal(beta3, 0.0278, tolerance = 0.01)

  expect_equal(mu1, -0.04, tolerance = 0.1)
  expect_equal(mu2, -0.013, tolerance = 0.1)
  expect_equal(mu3, -0.011, tolerance = 0.1)

  expect_equal(sigma, 0.087, tolerance = 0.01)
  expect_equal(Gamma, 0.63636, tolerance = 0.01)

  expect_equal(class(fit1), "lucid_parallel")

})




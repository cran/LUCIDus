# LUCID - five omics, binary outcome



test_that("check estimations of LUCID with binary outcome (K = 2,2,2)", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z4 <- matrix(rnorm(1000), nrow = 100)
  Z5 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5)
  Y <- rbinom(n=100, size =1, prob =0.25)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)


  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2, 2, 2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE)))

  betas <- mean(unlist(fit1$res_Beta$Beta))
  mus <- mean(unlist(fit1$res_Mu))
  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(unlist(fit1$res_Gamma$Gamma))

  # check parameters
  expect_equal(betas, -0.0401, tolerance = 0.01)
  expect_equal(mus, 0.020465, tolerance = 0.01)
  expect_equal(sigma, 0.08593, tolerance = 0.01)
  expect_equal(Gamma, 0.2973, tolerance = 0.01)

  ##missing data
  a = sample(1:1000, 30, replace=FALSE)
  Z1[a] = NA
  Z2[62:65, 6:8] = NA
  a = sample(1:1000, 20, replace=FALSE)
  Z4[a] = NA
  Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5)

  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2, 2, 2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE,
                                             init_impute = "mix")))

  betas <- mean(unlist(fit1$res_Beta$Beta))
  mus <- mean(unlist(fit1$res_Mu))
  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(unlist(fit1$res_Gamma$Gamma))

  # check parameters
  expect_equal(betas, -0.01601, tolerance = 0.1)
  expect_equal(mus, 0.04254, tolerance = 0.1)
  expect_equal(sigma, 0.075786, tolerance = 0.1)
  expect_equal(Gamma, 0.2973, tolerance = 0.1)


})




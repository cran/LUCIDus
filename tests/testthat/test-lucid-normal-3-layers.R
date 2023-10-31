# LUCID - three omics, normal outcome

test_that("check estimations of LUCID with normal outcome (K = 2,2,2)", {

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
                                                  init_omic.data.model  = "VVV",
                                                  seed = i,
                                                  init_impute = "mix",
                                                  init_par = "mclust",
                                                  useY = TRUE)))
  
  
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2), CoG = CoG, CoY = CoY,
                                                  lucid_model = "parallel",
                                                  family = "normal",
                                                  init_omic.data.model  = "VVV",
                                                  seed = i,
                                                  init_impute = "mix",
                                                  init_par = "mclust",
                                                  verbose = TRUE,
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
  expect_equal(beta1, 0.100, tolerance = 0.01)
  expect_equal(beta2, -0.236, tolerance = 0.01)
  expect_equal(beta3, -0.0256, tolerance = 0.01)

  expect_equal(mu1, -0.042, tolerance = 0.01)
  expect_equal(mu2, 0.1119, tolerance = 0.01)
  expect_equal(mu3, -0.01587, tolerance = 0.01)

  expect_equal(sigma, 0.07487, tolerance = 0.01)
  expect_equal(Gamma, 0.6765, tolerance = 0.01)

  expect_equal(class(fit1), "lucid_parallel")

})



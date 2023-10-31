# LUCID - LUCID in Serial, binary outcome



test_that("check estimations of LUCID in Serial with binary outcome (K = 2,2,2,2,2)", {
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

  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2, 2, 2),
                                             lucid_model = "serial",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             CoG = CoG, CoY = CoY,
                                             seed = i,
                                             useY = TRUE)))

  betas <- mean(unlist(fit1$res_Beta))
  mus <- mean(unlist(fit1$res_Mu))
  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(unlist(fit1$res_Gamma))

  # check parameters
  expect_equal(betas, -0.0672, tolerance = 0.01)
  expect_equal(mus, 0.0189, tolerance = 0.01)
  expect_equal(sigma, 0.086, tolerance = 0.01)
  expect_equal(Gamma, -0.164, tolerance = 0.01)

  expect_equal(class(fit1), "lucid_serial")


  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), Z4 = Z4, Z5 = Z5)
  invisible(capture.output(fit2 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2,list(2,2),2,2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "serial",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE)))

  betas <- mean(unlist(fit2$res_Beta))
  mus <- mean(unlist(fit2$res_Mu))
  sigma <- mean(unlist(fit2$res_Sigma))
  Gamma <- mean(unlist(fit2$res_Gamma))

  # check parameters
  expect_equal(betas, -0.0672, tolerance = 0.01)
  expect_equal(mus, 0.01795, tolerance = 0.01)
  expect_equal(sigma, 0.0863, tolerance = 0.01)
  expect_equal(Gamma, -0.164, tolerance = 0.01)

  expect_equal(class(fit2), "lucid_serial")


  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), list(Z4 = Z4, Z5 = Z5))
  invisible(capture.output(fit3 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(3,list(2,2),list(2,2)),
                                                  lucid_model = "serial",
                                                  family = "binary",
                                                  seed = i,
                                                  init_omic.data.model = "VVV",
                                                  CoG = CoG, CoY = CoY,
                                                  useY = TRUE)))

  betas <- mean(unlist(fit3$res_Beta))
  mus <- mean(unlist(fit3$res_Mu))
  sigma <- mean(unlist(fit3$res_Sigma))
  Gamma <- mean(unlist(fit3$res_Gamma$Gamma$mu))

  # check parameters
  expect_equal(betas, -7.73, tolerance = 0.01)
  expect_equal(mus, 0.01496, tolerance = 0.01)
  expect_equal(sigma, 0.0877, tolerance = 0.01)
  expect_equal(Gamma, 0.25, tolerance = 0.01)

  expect_equal(class(fit3), "lucid_serial")

})




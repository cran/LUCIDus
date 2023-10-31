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

  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2, 2, 2),
                                             lucid_model = "parallel",
                                             family = "binary",

                                             seed = i,
                                             useY = TRUE)))
  betas <- mean(unlist(fit1$res_Beta$Beta))
  mus <- mean(unlist(fit1$res_Mu))
  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(unlist(fit1$res_Gamma$Gamma))

  # check parameters
  expect_equal(betas, 0.03, tolerance = 0.05)


  expect_equal(mus, 0.02005, tolerance = 0.05)

  expect_equal(sigma, 0.09, tolerance = 0.05)
  expect_equal(Gamma, 0.297, tolerance = 0.05)
  expect_equal(class(fit1), "lucid_parallel")

  invisible(capture.output(fit2 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2, 2, 2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "binary",

                                             seed = i,
                                             useY = TRUE)))

  betas <- mean(unlist(fit2$res_Beta$Beta))
  mus <- mean(unlist(fit2$res_Mu))
  sigma <- mean(unlist(fit2$res_Sigma))
  Gamma <- mean(unlist(fit2$res_Gamma$Gamma))

  # check parameters
  expect_equal(betas, -0.0437, tolerance = 0.05)
  expect_equal(mus, 0.02005, tolerance = 0.05)
  expect_equal(sigma, 0.09, tolerance = 0.05)
  expect_equal(Gamma, 0.2973, tolerance = 0.05)
})




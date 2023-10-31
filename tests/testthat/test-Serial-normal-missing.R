# LUCID - LUCID in Serial, normal outcome

test_that("check estimations of LUCID in Serial with normal outcome (K = 2,list(2,2),3,2)", {
  # run LUCID model
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z4 <- matrix(rnorm(1000), nrow = 100)
  Z5 <- matrix(rnorm(1000), nrow = 100)

  ##missing data
  a = sample(1:1000, 30, replace=FALSE)
  Z1[a] = NA
  Z2[62:65, 6:8] = NA
  a = sample(1:1000, 30, replace=FALSE)
  Z4[a] = NA

  Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5)
  Y <- rnorm(100)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2,2,2,2,2),
                                             lucid_model = "serial",
                                             family = "normal",
                                             seed = i,
                                             CoG = CoG, CoY = CoY,
                                             useY = TRUE)))

  betas <- mean(unlist(fit1$res_Beta))
  mus <- mean(unlist(fit1$res_Mu))
  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(unlist(fit1$res_Gamma))



  expect_equal(class(fit1), "lucid_serial")

  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), Z4 = Z4, Z5 = Z5)
  invisible(capture.output(fit2 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2,list(2,2),2,2),
                                                  lucid_model = "serial",
                                                  family = "normal",
                                                  seed = i,
                                                  CoG = CoG, CoY = CoY,
                                                  useY = TRUE)))

  betas <- mean(unlist(fit2$res_Beta))
  mus <- mean(unlist(fit2$res_Mu))
  sigma <- mean(unlist(fit2$res_Sigma))
  Gamma <- mean(unlist(fit2$res_Gamma))



  expect_equal(class(fit2), "lucid_serial")


  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), list(Z4 = Z4, Z5 = Z5))
  invisible(capture.output(fit3 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(3,list(2,2),list(2,2)),
                                                  lucid_model = "serial",
                                                  family = "normal",
                                                  seed = i,
                                                  CoG = CoG, CoY = CoY,
                                                  useY = TRUE)))

  betas <- mean(unlist(fit3$res_Beta))
  mus <- mean(unlist(fit3$res_Mu))
  sigma <- mean(unlist(fit3$res_Sigma))
  Gamma <- mean(unlist(fit3$res_Gamma$Gamma$mu))



  expect_equal(class(fit3), "lucid_serial")


})

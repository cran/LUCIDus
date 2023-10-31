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

  sum_fit1 = summary_lucid(fit1)
  print.sumlucid(sum_fit1)


  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), Z4 = Z4, Z5 = Z5)
  invisible(capture.output(fit2 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2,list(2,2),2,2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "serial",
                                             family = "binary",

                                             seed = i,
                                             useY = TRUE)))

  sum_fit2 = summary_lucid(fit2)
  print.sumlucid(sum_fit2)


  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), list(Z4 = Z4, Z5 = Z5))
  invisible(capture.output(fit3 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(3,list(2,2),list(2,2)),
                                                  lucid_model = "serial",
                                                  family = "binary",
                                                  seed = i,
                                                  CoG = CoG, CoY = CoY,
                                                  useY = TRUE)))
  sum_fit3 = summary_lucid(fit3)
  print.sumlucid(sum_fit3)

})




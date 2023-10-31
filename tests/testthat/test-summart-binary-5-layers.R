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
  Z_e = matrix(rnorm(1000), nrow = 100)
  Y <- rbinom(n=100, size =1, prob =0.25)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)

  ## early
  invisible(capture.output(fit1 <- lucid(G = G, Z = Z_e, Y = Y, K = 2:4,
                                                  CoG = CoG, CoY = CoY,
                                                  lucid_model = "early",
                                                  family = "binary",init_omic.data.model = "VVV",
                                                  seed = i,
                                                  useY = TRUE)))

  sum_fit1 = summary_lucid(fit1)
  print.sumlucid(sum_fit1)



  invisible(capture.output(fit2 <- lucid(G = G, Z = Z, Y = Y, K = list(3, 2:4, 3, 2:3, 2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "binary",
                                         init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE)))
  sum_fit2 = summary_lucid(fit2)
  print.sumlucid(sum_fit2)


  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), Z4 = Z4, Z5 = Z5)
  invisible(capture.output(fit3 <- lucid(G = G, Z = Z, Y = Y, K = list(2:3,list(2,2:4),2:3,2),
                                                  CoG = CoG, CoY = CoY,
                                                  lucid_model = "serial",
                                                  family = "binary",
                                         init_omic.data.model = "VVV",
                                                  seed = i,
                                                  useY = TRUE)))

  sum_fit3 = summary_lucid(fit3)
  print.sumlucid(sum_fit3)

})




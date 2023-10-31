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


  Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5)
  Y <- rnorm(100)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- tune_lucid(G = G, Z = Z, Y = Y, K = list(2,2:4,2,2:3,2),
                                             lucid_model = "serial",
                                             family = "normal",
                                             seed = i,
                                             CoG = CoG, CoY = CoY,
                                             useY = TRUE)))

  invisible(capture.output(best1 <- lucid(G = G, Z = Z, Y = Y, K = list(2,2:4,2,2:3,2),
                                              lucid_model = "serial",
                                              family = "normal",
                                              seed = i,
                                              CoG = CoG, CoY = CoY,
                                              useY = TRUE)))

  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), Z4 = Z4, Z5 = Z5)
  invisible(capture.output(fit2 <- tune_lucid(G = G, Z = Z, Y = Y, K = list(2,list(2:3,2),2:4,2),
                                                  lucid_model = "serial",
                                                  family = "normal",
                                                  seed = i,
                                                  CoG = CoG, CoY = CoY,
                                                  useY = TRUE)))

  invisible(capture.output(best2 <- lucid(G = G, Z = Z, Y = Y, K = list(2,list(2:3,2),2:4,2),
                                              lucid_model = "serial",
                                              family = "normal",
                                              seed = i,
                                              CoG = CoG, CoY = CoY,
                                              useY = TRUE)))

  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), list(Z4 = Z4, Z5 = Z5))
  invisible(capture.output(fit3 <- tune_lucid(G = G, Z = Z, Y = Y, K = list(3:4,list(2,2),list(2:5,2)),
                                                  lucid_model = "serial",
                                                  family = "normal",
                                                  seed = i,
                                                  CoG = CoG, CoY = CoY,
                                                  useY = TRUE)))

  invisible(capture.output(best3 <- lucid(G = G, Z = Z, Y = Y, K = list(3:4,list(2,2),list(2:5,2)),
                                              lucid_model = "serial",
                                              family = "normal",
                                              seed = i,
                                              CoG = CoG, CoY = CoY,
                                              useY = FALSE,
                                              verbose_tune = TRUE)))


})

# LUCID - five omics, normal outcome

test_that("check the summary function of LUCID with normal outcome (K = 2,2,2,2,2)", {
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
  Z_e = matrix(rnorm(1000), nrow = 100)
  Y <- rnorm(100)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))

  ## early
  invisible(capture.output(fit1 <- lucid(G = G, Z = Z_e, Y = Y, K = 2:5,
                                                  CoG = CoG, CoY = CoY,
                                                  lucid_model = "early",
                                                  family = "normal",
                                                  seed = i,
                                                  useY = TRUE)))

  sum_fit1 = summary_lucid(fit1)
  print.sumlucid(sum_fit1)

  ## parallel
  invisible(capture.output(fit2 <- lucid(G = G, Z = Z, Y = Y, K = list(3, 2:3, 3:4, 2:3, 2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "normal",
                                             seed = i,
                                             useY = TRUE)))

  sum_fit2 = summary_lucid(fit2)
  print.sumlucid(sum_fit2)

  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), Z4 = Z4, Z5 = Z5)
  invisible(capture.output(fit3 <- lucid(G = G, Z = Z, Y = Y, K = list(2,list(2:3,2),2:4,2),
                                              lucid_model = "serial",
                                              family = "normal",
                                              seed = i,
                                              CoG = CoG, CoY = CoY,
                                              useY = TRUE)))
  sum_fit3 = summary_lucid(fit3)
  print.sumlucid(sum_fit3)

})

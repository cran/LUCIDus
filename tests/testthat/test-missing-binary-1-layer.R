# LUCID - 1 omics, binary outcome



test_that("check estimations of LUCID with binary outcome (K = 2)", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z <- matrix(rnorm(1000),nrow = 100)
  Y <- rbinom(n=100, size =1, prob =0.25)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)


  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = 2,
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "early",
                                             family = "binary",

                                             seed = i,
                                             useY = TRUE)))



})




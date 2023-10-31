# LUCID - three omics, binary outcome



test_that("check estimations of LUCID with binary outcome (K = 2,2,2)", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2)
  Y <- rbinom(n=100, size =1, prob =0.25)




  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(tune_list1 <- tune_lucid(G = G, Z = Z1, Y = Y, K = 2:5,CoG = CoG, CoY = CoY,
                                             lucid_model = "early",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE)))

  invisible(capture.output(best1 <- lucid(G = G, Z = Z1, Y = Y, K = 2:5,CoG = CoG, CoY = CoY,
                                                    lucid_model = "early",
                                                    family = "binary",
                                          init_omic.data.model = "VVV",
                                                    seed = i,
                                                    useY = TRUE)))


  invisible(capture.output(tune_list2 <- tune_lucid(G = G, Z = Z, Y = Y, K = list(2:4,2),
                                                    CoG = CoG, CoY = CoY,
                                             lucid_model = "parallel",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE)))

  invisible(capture.output(tune_list2 <- tune_lucid(G = G, Z = Z, Y = Y, K = list(4,2),
                                                    CoG = CoG, CoY = CoY,
                                                    lucid_model = "parallel",
                                                    family = "binary",
                                                    init_omic.data.model = "VVV",
                                                    seed = i,
                                                    useY = TRUE)))

  invisible(capture.output(best2 <- lucid(G = G, Z = Z, Y = Y, K = list(2:4,2),
                                                    CoG = CoG, CoY = CoY,
                                                    lucid_model = "parallel",
                                                    family = "binary",
                                          init_omic.data.model = "VVV",
                                                    seed = i,
                                                    useY = TRUE)))
  invisible(capture.output(best3 <- lucid(G = G, Z = Z, Y = Y, K = list(4,2),
                                          CoG = CoG, CoY = CoY,
                                          lucid_model = "parallel",
                                          family = "binary",
                                          init_omic.data.model = "VVV",
                                          seed = i,
                                          useY = TRUE)))
})




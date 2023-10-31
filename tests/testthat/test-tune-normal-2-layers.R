# LUCID - three omics, normal outcome

test_that("check estimations of LUCID with normal outcome (K = 2,2,2)", {

  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  Y <- rnorm(100)
  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))

  invisible(capture.output(tune_list1 <- tune_lucid(G = G, Z = Z2, Y = Y, K = 2:4,
                                              CoG = CoG, CoY = CoY,
                                              lucid_model = "early",
                                              family = "normal",
                                              init_omic.data.model = "VVV",
                                              seed = i,
                                              init_impute = "mix",
                                              init_par = "mclust",
                                              useY = TRUE)))
  invisible(capture.output(best1 <- lucid(G = G, Z = Z2, Y = Y, K = 2:4,
                                                    CoG = CoG, CoY = CoY,
                                                    lucid_model = "early",
                                                    family = "normal",
                                                    init_omic.data.model = "VVV",
                                                    seed = i,
                                                    init_impute = "mix",
                                                    init_par = "mclust",
                                                    useY = TRUE)))

  invisible(capture.output(tune_list2 <- tune_lucid(G = G, Z = Z, Y = Y, K = list(2:4,2),
                                                  CoG = CoG, CoY = CoY,
                                                  lucid_model = "parallel",
                                                  family = "normal",
                                                  init_omic.data.model = "VVV",
                                                  seed = i,
                                                  init_impute = "mix",
                                                  init_par = "mclust",
                                                  useY = TRUE)))
  invisible(capture.output(best2 <- lucid(G = G, Z = Z, Y = Y, K = list(2:4,2),
                                                    CoG = CoG, CoY = CoY,
                                                    lucid_model = "serial",
                                                    family = "normal",
                                                    init_omic.data.model = "VVV",
                                                    seed = i,
                                                    init_impute = "mix",
                                                    init_par = "mclust",
                                                    useY = TRUE)))

})



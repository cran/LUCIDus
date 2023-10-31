# LUCID - two omics, normal outcome

test_that("check estimations of LUCID with normal outcome (K = 2,2)", {
  # run LUCID model
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z <- matrix(rnorm(1000),nrow = 100)
  a = sample(1:1000, 10, replace=FALSE)
  Z[a] = NA
  Z[66, 7:8] = NA
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  Y <- rnorm(100)
  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = 2,CoG = CoG, CoY = CoY,
                                             lucid_model = "early",
                                             family = "normal",

                                             seed = i,
                                             useY = TRUE)))




})



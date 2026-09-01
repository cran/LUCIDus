# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# LUCID - three omics, normal outcome

test_that("check estimations of LUCID with normal outcome (K = 2,2,2)", {
  skip_if_not_installed("mix")

  # run LUCID model
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2, Z2 = Z3)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  Y <- rnorm(100)


  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2), CoG = CoG, CoY = CoY,
                                                  lucid_model = "parallel",
                                                  family = "normal",
                                                  init_omic.data.model  = "VVV",
                                                  seed = i,
                                                  init_impute = "mix",
                                                  init_par = "mclust",
                                                  useY = TRUE)))
  
  
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2), CoG = CoG, CoY = CoY,
                                                  lucid_model = "parallel",
                                                  family = "normal",
                                                  init_omic.data.model  = "VVV",
                                                  seed = i,
                                                  init_impute = "mix",
                                                  init_par = "mclust",
                                                  verbose = TRUE,
                                                  useY = TRUE)))
  
  betas <- fit1$res_Beta$Beta
  beta1 <- mean(unlist(betas[1]))
  beta2 <- mean(unlist(betas[2]))
  beta3 <- mean(unlist(betas[3]))

  mus <- fit1$res_Mu
  mu1 <- mean(unlist(mus[1]))
  mu2 <- mean(unlist(mus[2]))
  mu3 <- mean(unlist(mus[3]))

  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(parallel_delta_coef(fit1$res_Gamma$Gamma))

  # check parameters via robust invariants
  expect_true(all(is.finite(c(beta1, beta2, beta3))))
  expect_true(all(is.finite(c(mu1, mu2, mu3))))
  expect_true(is.finite(sigma))
  expect_true(is.finite(Gamma))
  expect_gt(sigma, 0)
  expect_true(all(abs(c(beta1, beta2, beta3)) < 2))
  expect_true(all(abs(c(mu1, mu2, mu3)) < 2))

  expect_equal(class(fit1), "lucid_parallel")

})

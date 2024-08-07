# LUCID - one omics, normal outcome


test_that("check estimations of LUCID with normal outcome (K = 2)", {
  # run LUCID model
  G <- sim_data$G[1:200, ]
  Z <- sim_data$Z[1:200, ]
  Y_normal <- sim_data$Y_normal[1:200, ]
  cov <- sim_data$Covariate[1:200, ]
  X <- sim_data$X[1:200]
  # i <- sample(1:2000, 1)
  i <- 1008
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- lucid(G = G,
                                         Z = Z,
                                         Y = Y_normal,
                                         CoY = cov,
                                         lucid_model = "early",
                                         family = "normal",
                                         K = 2,
                                         seed = i)))
  pars <- fit1
  beta_causal <- mean(pars$res_Beta[2, 2:5])
  beta_non <- mean(pars$res_Beta[2, 6:10])
  mu_causal <- mean(abs(pars$res_Mu[1, 1:5] - pars$res_Mu[2, 1:5]))
  mu_non <- mean(abs(pars$res_Mu[1, 6:10] - pars$res_Mu[2, 6:10]))
  gamma_causal <- as.numeric(abs(pars$res_Gamma$beta[1] - pars$res_Gamma$beta[2]))
  gamma_non <- as.numeric(mean(pars$res_Gamma$beta[3:4]))
  sigma <- mean(unlist(fit1$res_Sigma))

  # check parameters
  expect_equal(beta_causal, log(2), tolerance = 0.2)
  expect_equal(beta_non, 0, tolerance = 0.1)
  expect_equal(mu_causal, 2, tolerance = 0.1)
  expect_equal(mu_non, 0, tolerance = 0.1)
  expect_equal(gamma_causal, 1, tolerance = 0.05)
  expect_equal(gamma_non, 0, tolerance = 0.05)
  expect_equal(sigma, 0.1048542, tolerance = 0.05)

  # check summary_lucid
  sum_fit1 <- summary(fit1)
  expect_equal(class(fit1), "early_lucid")
  expect_equal(class(sum_fit1), "sumlucid_early")
})


test_that("check variable selection on G", {
  # run LUCID model
  G <- sim_data$G[1:200, ]
  Z <- sim_data$Z[1:200, ]
  Y_normal <- sim_data$Y_normal[1:200, ]
  cov <- sim_data$Covariate[1:200, ]
  X <- sim_data$X[1:200]
  # i <- sample(1:2000, 1)
  i <- 1008
  # cat(paste("test2 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- lucid(G = G,
                                         Z = Z,
                                         Y = Y_normal,
                                         CoY = cov,
                                         lucid_model = "early",
                                         family = "normal",
                                         K = 2,
                                         seed = i,

                                         Rho_G = 0.05)))

  # check parameters
  expect_equal(class(fit1$select$selectG), "logical")
  expect_equal(as.vector(fit1$select$selectG),
               rep(TRUE, 4))
})


test_that("check variable selection on Z", {
  # run LUCID model
  G <- sim_data$G[1:200, ]
  Z <- sim_data$Z[1:200, ]
  Y_normal <- sim_data$Y_normal[1:200, ]
  cov <- sim_data$Covariate[1:200, ]
  X <- sim_data$X[1:200]
  # i <- sample(1:2000, 1)
  i <- 1008
  # cat(paste("test3 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- lucid(G = G,
                                         Z = Z,
                                         Y = Y_normal,
                                         CoY = cov,
                                         lucid_model = "early",
                                         family = "normal",
                                         K = 2,
                                         seed = i,

                                         init_par = "random",
                                         Rho_Z_Mu = 13,
                                         Rho_Z_Cov = 0.05)))

  # check parameters
  expect_equal(class(fit1$select$selectG), "logical")
})


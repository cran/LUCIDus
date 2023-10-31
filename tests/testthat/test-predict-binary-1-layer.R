# LUCID - 1 omics, binary outcome

test_that("check estimations of LUCID with binary outcome (K = 2)", {


  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  i <- 1008

  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = 2,
                                             lucid_model = "early",
                                             family = "binary",
                                             CoY = cov,
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE)))


  pars <- fit1
  beta_causal <- mean(pars$res_Beta[2, 2:5])
  beta_non <- mean(pars$res_Beta[2, 6:10])
  mu_causal <- mean(abs(pars$res_Mu[1, 1:5] - pars$res_Mu[2, 1:5]))
  mu_non <- mean(abs(pars$res_Mu[1, 6:10] - pars$res_Mu[2, 6:10]))
  gamma <- as.numeric(pars$res_Gamma$beta)

  # check parameters
  expect_equal(beta_causal, log(2), tolerance = 0.2)
  expect_equal(beta_non, 0, tolerance = 0.1)
  expect_equal(mu_causal, 2, tolerance = 0.1)
  expect_equal(mu_non, 0, tolerance = 0.1)
  expect_equal(gamma, c(-0.5, 0.9, 0.8, -0.8), tolerance = 0.2)

  # check summary_lucid
  sum_fit1 <- summary_lucid(fit1)
  expect_equal(class(fit1), "early_lucid")
  expect_equal(class(sum_fit1), "sumlucid_early")
})

test_that("check variable selection on G", {
  # run LUCID model
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  # cat(paste("test2 - seed =", i, "\n"))
  i <- 1008
  invisible(capture.output(fit1 <- estimate_lucid(G = G,
                                         Z = Z,
                                         Y = Y_binary,
                                         CoY = cov,
                                         family = "binary",
                                         lucid_model = "early",
                                         K = 2,
                                         seed = i,
                                         useY = TRUE

                                          )))

  # make prediction
  pred1 <- predict_lucid(model = fit1,
                         lucid_model = "early",
                         G = G,
                         Z = Z,
                         Y = Y_binary,
                         CoY = cov)


  # compare prediction of X
  expect_equal(fit1$inclusion.p, pred1$inclusion.p, tolerance = 0.05)
  expect_equal(class(pred1$pred.x), "integer")

  pred2 <- predict_lucid(model = fit1,
                         lucid_model = "early",
                         G = G,
                         Z = Z,
                         Y = NULL,
                         CoY = cov)

  expect_equal(class(pred2$pred.x), "integer")

})

# LUCID - 1 omics, binary outcome

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
  expect_equal(class(pred1$pred.x), "numeric")

  pred2 <- predict_lucid(model = fit1,
                         lucid_model = "early",
                         G = G,
                         Z = Z,
                         Y = NULL,
                         CoY = cov)

  expect_equal(class(pred2$pred.x), "numeric")
  

  counter_G_v1_1 <- G
  counter_G_v1_1[,1] = 1
  counter_G_v1_0 <- G
  counter_G_v1_0[,1] = 0
  g_comp_1 <- predict_lucid(model = fit1,
                         lucid_model = "early",
                         G = counter_G_v1_1,
                         Z = Z,
                         Y = NULL,
                         g_computation = TRUE,
                         CoY = cov)
  g_comp_0 <- predict_lucid(model = fit1,
                            lucid_model = "early",
                            G = counter_G_v1_0,
                            Z = Z,
                            Y = NULL,
                            g_computation = TRUE,
                            CoY = cov)
  
})

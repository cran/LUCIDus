# test prediction function for LUCID


test_that("check prediction of lucid", {
  # run LUCID model
  G <- sim_data$G[1:200, ]
  Z <- sim_data$Z[1:200, ]
  Y_normal <- sim_data$Y_normal[1:200, ]
  cov <- sim_data$Covariate[1:200, ]
  X <- sim_data$X[1:200]
  i <- 1008
  invisible(capture.output(fit1 <- estimate_lucid(G = G,
                                             Z = Z,
                                             Y = Y_normal,
                                             CoY = cov,
                                             lucid_model = "early",
                                             init_omic.data.model = "VVV",
                                             family = "normal",
                                             K = 2,
                                             seed = i)))

  # make prediction
  pred1 <- predict_lucid(model = fit1,
                         lucid_model = "early",
                         G = G,
                         Z = Z,
                         Y = Y_normal,
                         CoY = cov)

  dat <- as.data.frame(cbind(Y_normal, fit1$inclusion.p[, -1], cov))
  fit_lm <- lm(Y_normal ~., data = dat)
  pred2 <- as.vector(predict(fit_lm))
  # compare prediction of X
  expect_equal(fit1$inclusion.p, pred1$inclusion.p)
  expect_equal(class(pred1$pred.x), "integer")
  expect_equal(pred1$pred.y, pred2)

  pred2 <- predict_lucid(model = fit1,
                         lucid_model = "early",
                         G = G,
                         Z = Z,

                         Y = NULL,
                         CoY = cov)

  expect_equal(class(pred2$pred.x), "integer")
})

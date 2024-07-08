# LUCID - five omics, normal outcome

test_that("check predictions of LUCID with normal outcome (K = 2,2,2,2,2)", {
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
  Y <- rnorm(100)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2, 2, 2, 2, 2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "serial",
                                             init_omic.data.model = "VVV",
                                             family = "normal",
                                             seed = i,
                                             useY = TRUE)))
  #use training data
  pred1 <- predict_lucid(model = fit1,
                         lucid_model = "serial",
                         G = G,
                         Z = Z,
                         Y = Y,
                         CoG = CoG, CoY = CoY)
  expect_equal(fit1$inclusion.p, pred1$inclusion.p, tolerance = 0.05)
  expect_equal(class(pred1$pred.x), "list")


  set.seed(i+1000)
  n_G <- matrix(rnorm(500), nrow = 100)
  n_Z1 <- matrix(rnorm(1000),nrow = 100)
  n_Z2 <- matrix(rnorm(1000), nrow = 100)
  n_Z3 <- matrix(rnorm(1000), nrow = 100)
  n_Z4 <- matrix(rnorm(1000), nrow = 100)
  n_Z5 <- matrix(rnorm(1000), nrow = 100)
  n_Z <- list(Z1 = n_Z1, Z2 = n_Z2, Z3 = n_Z3, Z4 = n_Z4, Z5 = n_Z5)
  n_Y <- rnorm(100)
  n_CoY <- matrix(rnorm(200), nrow = 100)
  n_CoG <- matrix(rnorm(200), nrow = 100)

  #use new data
  pred2 <- predict_lucid(model = fit1,
                         lucid_model = "serial",
                         G = n_G,
                         Z = n_Z,
                         Y = n_Y,
                         CoG = n_CoG, CoY = n_CoY)

  expect_equal(class(pred2$pred.x), "list")
  expect_equal(mean(pred2$pred.y), 0.1434522, tolerance = 0.05)
  expect_equal(mean(pred2$inclusion.p[[1]]), 0.5)

  pred3 <- predict_lucid(model = fit1,
                         lucid_model = "serial",
                         G = n_G,
                         Z = n_Z,
                         Y = NULL,
                         CoG = n_CoG, CoY = n_CoY)

  expect_equal(class(pred3$pred.x), "list")
  expect_equal(mean(pred3$pred.y), 0.1429235, tolerance = 0.05)
  expect_equal(mean(pred3$inclusion.p[[1]]), 0.5)
  
  G_gcomp_1 = G
  G_gcomp_1[,4] = 1
  G_gcomp_0 = G
  G_gcomp_0[,4] = 0
  
  gcomp1 <- predict_lucid(model = fit1,
                         lucid_model = "serial",
                         G = G_gcomp_1,
                         Z = Z,
                         Y = NULL,
                         g_computation = TRUE,
                         CoG = CoG, CoY = CoY)
  gcomp0 <- predict_lucid(model = fit1,
                          lucid_model = "serial",
                          G = G_gcomp_0,
                          Z = Z,
                          Y = NULL,
                          g_computation = TRUE,
                          CoG = CoG, CoY = CoY)
  
  Z <- list(list(Z1 = Z1, Z2 = Z2), Z3 = Z3, Z4 = Z4, Z5 = Z5)
  n_Z <- list(list(Z1 = n_Z1, Z2 = n_Z2), Z3 = n_Z3, Z4 = n_Z4, Z5 = n_Z5)
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(list(2, 2), 2, 2, 2),
                                                  CoG = CoG, CoY = CoY,
                                                  lucid_model = "serial",
                                                  family = "normal",
                                                  seed = i,
                                                  useY = TRUE)))

  #use training data
  pred1 <- predict_lucid(model = fit1,
                         lucid_model = "serial",
                         G = G,
                         Z = Z,
                         Y = Y,
                         CoG = CoG, CoY = CoY)
  expect_equal(fit1$inclusion.p, pred1$inclusion.p, tolerance = 0.05)
  expect_equal(class(pred1$pred.x), "list")


  #use new data
  pred2 <- predict_lucid(model = fit1,
                         lucid_model = "serial",
                         G = n_G,
                         Z = n_Z,
                         Y = n_Y,
                         CoG = n_CoG, CoY = n_CoY)

  expect_equal(class(pred2$pred.x), "list")
  expect_equal(mean(pred2$pred.y), 0.05, tolerance = 0.05)
  expect_equal(mean(pred2$inclusion.p[[2]]), 0.5)

  pred3 <- predict_lucid(model = fit1,
                         lucid_model = "serial",
                         G = n_G,
                         Z = n_Z,
                         Y = NULL,
                         CoG = n_CoG, CoY = n_CoY)

  expect_equal(class(pred3$pred.x), "list")
  expect_equal(mean(pred3$pred.y), 0.05, tolerance = 0.05)
  expect_equal(mean(pred3$inclusion.p[[2]]), 0.5)
  
  gcomp1 <- predict_lucid(model = fit1,
                          lucid_model = "serial",
                          G = G_gcomp_1,
                          Z = Z,
                          Y = NULL,
                          g_computation = TRUE,
                          CoG = CoG, CoY = CoY)
  gcomp0 <- predict_lucid(model = fit1,
                          lucid_model = "serial",
                          G = G_gcomp_0,
                          Z = Z,
                          Y = NULL,
                          g_computation = TRUE,
                          CoG = CoG, CoY = CoY)
  
})

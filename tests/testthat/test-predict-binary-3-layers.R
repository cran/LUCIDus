# LUCID - three omics, binary outcome



test_that("check prediction of LUCID with binary outcome (K = 2,2,2)", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3)
  Y <- rbinom(n=100, size =1, prob =0.65)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)

  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2),CoG = CoG, CoY = CoG,
                                             lucid_model = "parallel",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             useY = TRUE)))
  set.seed(i+1000)
  n_G <- matrix(rnorm(500), nrow = 100)
  n_Z1 <- matrix(rnorm(1000),nrow = 100)
  n_Z2 <- matrix(rnorm(1000), nrow = 100)
  n_Z3 <- matrix(rnorm(1000), nrow = 100)
  n_Z <- list(Z1 = n_Z1, Z2 = n_Z2, Z3 = n_Z3)
  n_Y <- rbinom(n=100, size =1, prob =0.25)
  n_CoY <- matrix(rnorm(200), nrow = 100)
  n_CoG <- matrix(rnorm(200), nrow = 100)

  #use training data
  pred1 <- predict_lucid(model = fit1,
                         lucid_model = "parallel",
                         G = G,
                         Z = Z,
                         Y = Y,
                         CoG = CoG, CoY = CoY, response = TRUE)

  expect_equal(fit1$inclusion.p, pred1$inclusion.p, tolerance = 0.05)
  expect_equal(class(pred1$pred.x), "list")
  expect_equal(max(pred1$pred.y), 1)
  expect_equal(mean(pred1$pred.y), 0.99)
  expect_equal(mean(pred1$inclusion.p[[1]]), 0.5)

  #use new data
  pred2 <- predict_lucid(model = fit1,
                         lucid_model = "parallel",
                         G = n_G,
                         Z = n_Z,
                         Y = n_Y,
                         CoG = n_CoG, CoY = n_CoY, response = TRUE)

  expect_equal(class(pred2$pred.x), "list")
  expect_equal(max(pred2$pred.y), 1)
  expect_equal(mean(pred2$pred.y), 1)
  expect_equal(mean(pred2$inclusion.p[[1]]), 0.5)

  #new data not using Y, and response = FALSE
  pred3 <- predict_lucid(model = fit1,
                         lucid_model = "parallel",
                         G = n_G,
                         Z = n_Z,
                         Y = NULL,
                         CoG = n_CoG, CoY = n_CoY, response = FALSE)

  expect_equal(class(pred3$pred.x), "list")
  expect_equal(max(pred3$pred.y), 0.9124988, tolerance = 0.05)
  expect_equal(mean(pred3$pred.y), 0.7277153, tolerance = 0.05)
  expect_equal(mean(pred3$inclusion.p[[1]]), 0.5)

})




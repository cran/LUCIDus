# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# LUCID - five omics, binary outcome



test_that("check estimations of LUCID with binary outcome (K = 2,2,2)", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z4 <- matrix(rnorm(1000), nrow = 100)
  Z5 <- matrix(rnorm(1000), nrow = 100)
  Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5)
  Y <- rbinom(n=100, size =1, prob =0.45)

  #dont use Cog Coy here

  invisible(capture.output(fit1 <- est_lucid(G = G, Z = Z, Y = Y, K = c(2, 2, 2, 2, 2),
                                             lucid_model = "parallel",
                                             family = "binary",

                                             seed = i,
                                             useY = TRUE)))
  set.seed(i+1000)
  n_G <- matrix(rnorm(500), nrow = 100)
  n_Z1 <- matrix(rnorm(1000),nrow = 100)
  n_Z2 <- matrix(rnorm(1000), nrow = 100)
  n_Z3 <- matrix(rnorm(1000), nrow = 100)
  n_Z4 <- matrix(rnorm(1000), nrow = 100)
  n_Z5 <- matrix(rnorm(1000), nrow = 100)
  n_Z <- list(Z1 = n_Z1, Z2 = n_Z2, Z3 = n_Z3, Z4 = n_Z4, Z5 = n_Z5)
  n_Y <- rbinom(n=100, size =1, prob =0.45)

  #use training data
  pred1 <- predict_lucid(model = fit1,
                         lucid_model = "parallel",
                         G = G,
                         Z = Z,
                         Y = Y,
                         response = TRUE)

  expect_equal(fit1$inclusion.p, pred1$inclusion.p, tolerance = 0.05)
  expect_true(is.list(pred1$pred.x))
  expect_equal(length(pred1$pred.x), length(fit1$K))
  expect_true(all(pred1$pred.y %in% c(0, 1)))
  expect_equal(length(pred1$pred.y), nrow(G))
  expect_true(all(sapply(seq_along(pred1$inclusion.p), function(j) {
    all(is.finite(pred1$inclusion.p[[j]])) &&
      all(rowSums(pred1$inclusion.p[[j]]) > 0.9999) &&
      all(rowSums(pred1$inclusion.p[[j]]) < 1.0001)
  })))

  #use new data
  pred2 <- predict_lucid(model = fit1,
                         lucid_model = "parallel",
                         G = n_G,
                         Z = n_Z,
                         Y = n_Y,
                         response = TRUE)

  expect_true(is.list(pred2$pred.x))
  expect_equal(length(pred2$pred.x), length(fit1$K))
  expect_true(all(pred2$pred.y %in% c(0, 1)))
  expect_equal(length(pred2$pred.y), nrow(n_G))
  expect_true(all(sapply(seq_along(pred2$inclusion.p), function(j) {
    all(is.finite(pred2$inclusion.p[[j]])) &&
      all(rowSums(pred2$inclusion.p[[j]]) > 0.9999) &&
      all(rowSums(pred2$inclusion.p[[j]]) < 1.0001)
  })))

  #new data not using Y, and response = FALSE
  pred3 <- predict_lucid(model = fit1,
                         lucid_model = "parallel",
                         G = n_G,
                         Z = n_Z,
                         Y = NULL,
                         response = FALSE)

  expect_true(is.list(pred3$pred.x))
  expect_equal(length(pred3$pred.x), length(fit1$K))
  expect_equal(length(pred3$pred.y), nrow(n_G))
  expect_true(all(is.finite(pred3$pred.y)))
  expect_true(all(pred3$pred.y >= 0 & pred3$pred.y <= 1))
  expect_true(all(sapply(seq_along(pred3$inclusion.p), function(j) {
    all(is.finite(pred3$inclusion.p[[j]])) &&
      all(rowSums(pred3$inclusion.p[[j]]) > 0.9999) &&
      all(rowSums(pred3$inclusion.p[[j]]) < 1.0001)
  })))

})



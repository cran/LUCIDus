test_that("Equation 17 reduces to the Gaussian conditional optimum for one cluster", {
  obs <- c(2, NA)
  mu <- matrix(c(0, 0), nrow = 2, ncol = 1)
  sigma <- array(matrix(c(1, 0.8, 0.8, 1), 2, 2), dim = c(2, 2, 1))

  got <- fill_data(obs, mu, sigma, p = 1, index = c(TRUE, FALSE),
                   lucid_model = "early")

  expect_equal(got[1], 2)
  expect_equal(got[2], 0.8 * 2, tolerance = 1e-12)
})

test_that("Equation 17 uses responsibility-density and precision weights", {
  obs <- c(0.7, -0.2, 0.4)
  index <- c(TRUE, FALSE, FALSE)
  mu_km <- rbind(c(-0.3, 0.5, -0.4), c(0.8, -0.6, 0.9))
  sigma <- array(0, c(3, 3, 2))
  sigma[, , 1] <- matrix(c(1.2, .3, -.2, .3, .8, .25, -.2, .25, 1.1), 3)
  sigma[, , 2] <- matrix(c(.9, -.15, .2, -.15, 1.4, .35, .2, .35, 1.0), 3)
  p <- c(.35, .65)
  A <- which(index)
  B <- which(!index)
  precision <- lapply(1:2, function(k) solve(sigma[, , k]))
  log_w <- log(p) + vapply(1:2, function(k) {
    mclust::dmvnorm(matrix(obs, nrow = 1), mu_km[k, ], sigma[, , k], log = TRUE)
  }, numeric(1))
  w <- exp(log_w - max(log_w)); w <- w / sum(w)
  lhs <- Reduce(`+`, lapply(1:2, function(k) w[k] * precision[[k]][B, B, drop = FALSE]))
  rhs <- Reduce(`+`, lapply(1:2, function(k) {
    omega <- precision[[k]]
    w[k] * (omega[B, A, drop = FALSE] %*% mu_km[k, A] +
              omega[B, B, drop = FALSE] %*% mu_km[k, B] -
              omega[B, A, drop = FALSE] %*% obs[A])
  }))
  expected <- solve(lhs, rhs)

  got <- fill_data(obs, mu_km, sigma, p, index, lucid_model = "parallel")
  expect_equal(got[B], as.numeric(expected), tolerance = 1e-10)
  expect_equal(got[A], obs[A])
})

test_that("I-step updates initialized sporadic cells using the original mask", {
  Z <- matrix(c(2, 0, -1, 0), nrow = 2, byrow = TRUE)
  index <- matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2, byrow = TRUE)
  mu <- matrix(c(0, 0), nrow = 1)
  sigma <- array(matrix(c(1, .75, .75, 1), 2), dim = c(2, 2, 1))
  p <- matrix(1, 2, 1)

  got <- Istep_Z(Z, p, mu, sigma, index, lucid_model = "parallel")
  expect_equal(got[, 2], .75 * Z[, 1], tolerance = 1e-12)
})

test_that("listwise-missing rows remain excluded rather than dynamically imputed", {
  Z <- matrix(c(1, 2, NA, NA), nrow = 2, byrow = TRUE)
  index <- !is.na(Z)
  mu <- matrix(c(0, 0), nrow = 1)
  sigma <- array(diag(2), dim = c(2, 2, 1))

  got_early <- Istep_Z(Z, matrix(1, 2, 1), mu, sigma, index,
                       lucid_model = "early")
  got_parallel <- Istep_Z(Z, matrix(1, 2, 1), mu, sigma, index,
                          lucid_model = "parallel")

  expect_true(all(is.na(got_early[2, ])))
  expect_true(all(is.na(got_parallel[2, ])))
  expect_equal(got_early[1, ], Z[1, ])
  expect_equal(got_parallel[1, ], Z[1, ])
})

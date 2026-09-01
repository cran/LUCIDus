skip_on_cran()

simulate_early_paper <- function(seed = 9201, N = 650) {
  set.seed(seed)
  G <- cbind(g1 = rnorm(N), g2 = rnorm(N))
  X <- rbinom(N, 1, plogis(-0.25 + G[, 1] - 0.45 * G[, 2])) + 1L
  Z <- cbind(
    z1 = rnorm(N, c(-1.5, 1.5)[X], 0.65),
    z2 = rnorm(N, c(-1, 1)[X], 0.70)
  )
  age <- rnorm(N)
  CoY <- matrix(age, ncol = 1, dimnames = list(NULL, "age"))
  list(
    G = G, X = X, Z = Z, CoY = CoY,
    Y_binary = rbinom(N, 1, plogis(c(-1.1, 0.9)[X] + 0.5 * age)),
    Y_normal = c(-0.8, 1.3)[X] + 0.6 * age + rnorm(N, sd = c(0.6, 1.1)[X])
  )
}

test_that("early LUCID recovers generating parameters for both outcome families", {
  d <- simulate_early_paper()
  fits <- lapply(c("binary", "normal"), function(family) {
    fit <- NULL
    invisible(capture.output(fit <- est_lucid(
      "early", d$G, d$Z, d[[paste0("Y_", family)]], CoY = d$CoY,
      K = 2, init_omic.data.model = "VVV", family = family,
      seed = 10, max_itr = 100, max_tot.itr = 300
    )))
    fit
  })
  binary_fit <- fits[[1]]
  normal_fit <- fits[[2]]

  expect_gt(mean(max.col(binary_fit$inclusion.p) == d$X), 0.97)
  expect_equal(unname(binary_fit$res_Beta[2, ]), c(-0.25, 1, -0.45), tolerance = 0.30)
  expect_equal(unname(binary_fit$res_Gamma$beta), c(-1.1, 2, 0.5), tolerance = 0.30)
  expect_equal(binary_fit$res_Mu[, 1], c(-1.5, 1.5), tolerance = 0.20)

  expect_gt(mean(max.col(normal_fit$inclusion.p) == d$X), 0.97)
  expect_equal(unname(normal_fit$res_Gamma$beta), c(-0.8, 2.1, 0.6), tolerance = 0.20)
  expect_equal(unname(normal_fit$res_Gamma$sigma), c(0.6, 1.1), tolerance = 0.15)
})

simulate_parallel_paper <- function(seed = 9202, N = 700) {
  set.seed(seed)
  G <- cbind(g1 = rnorm(N), g2 = rnorm(N))
  X1 <- rbinom(N, 1, plogis(-0.2 + 0.9 * G[, 1])) + 1L
  X2 <- rbinom(N, 1, plogis(0.25 - 0.8 * G[, 2])) + 1L
  Z <- list(
    cbind(a1 = rnorm(N, c(-1.3, 1.3)[X1], .6),
          a2 = rnorm(N, c(-.8, .8)[X1], .65)),
    cbind(b1 = rnorm(N, c(-1.4, 1.4)[X2], .6),
          b2 = rnorm(N, c(-.7, .7)[X2], .7))
  )
  age <- rnorm(N)
  CoY <- matrix(age, ncol = 1, dimnames = list(NULL, "age"))
  list(
    G = G, X1 = X1, X2 = X2, Z = Z, CoY = CoY,
    Y_normal = -0.5 + 0.9 * (X1 == 2) + 0.6 * (X2 == 2) +
      0.4 * age + rnorm(N, sd = 0.75),
    Y_binary = rbinom(N, 1, plogis(-1 + 0.8 * (X1 == 2) +
                                      0.55 * (X2 == 2) + 0.35 * age))
  )
}

test_that("parallel LUCID recovers additive effects for normal and binary outcomes", {
  d <- simulate_parallel_paper()
  fit_model <- function(family) {
    fit <- NULL
    invisible(capture.output(fit <- est_lucid(
      "parallel", d$G, d$Z, d[[paste0("Y_", family)]], CoY = d$CoY,
      K = c(2, 2), init_omic.data.model = "VVV", family = family,
      seed = 12, max_itr = 100, max_tot.itr = 300
    )))
    fit
  }
  normal_fit <- fit_model("normal")
  binary_fit <- fit_model("binary")

  expect_gt(mean(max.col(normal_fit$inclusion.p[[1]]) == d$X1), 0.97)
  expect_gt(mean(max.col(normal_fit$inclusion.p[[2]]) == d$X2), 0.97)
  expect_equal(unname(parallel_delta_coef(normal_fit$res_Gamma$Gamma)),
               c(-0.5, 0.9, 0.6, 0.4), tolerance = 0.20)
  expect_equal(normal_fit$res_Gamma$Gamma$sd, 0.75, tolerance = 0.10)
  expect_equal(unname(parallel_delta_coef(binary_fit$res_Gamma$Gamma)),
               c(-1, 0.8, 0.55, 0.35), tolerance = 0.30)
})

test_that("serial LUCID recovers sequential clusters and final outcome effect", {
  set.seed(9205)
  N <- 600
  G <- cbind(g1 = rnorm(N), g2 = rnorm(N))
  X1 <- rbinom(N, 1, plogis(-0.2 + 0.9 * G[, 1])) + 1L
  X2 <- rbinom(N, 1, plogis(-0.5 + 1.8 * (X1 == 2))) + 1L
  Z <- list(
    cbind(rnorm(N, c(-1.4, 1.4)[X1], .6), rnorm(N, c(-1, 1)[X1], .65)),
    cbind(rnorm(N, c(-1.5, 1.5)[X2], .6), rnorm(N, c(-.9, .9)[X2], .65))
  )
  Y <- c(-0.8, 1.2)[X2] + rnorm(N, sd = 0.7)
  fit <- NULL
  invisible(capture.output(fit <- estimate_lucid(
    "serial", G, Z, Y, K = list(2, 2), init_omic.data.model = "VVV",
    family = "normal", seed = 15, max_itr = 100, max_tot.itr = 300
  )))
  acc1 <- mean(max.col(fit$submodel[[1]]$inclusion.p) == X1)

  expect_gt(max(acc1, 1 - acc1), 0.97)
  expect_gt(mean(max.col(fit$submodel[[2]]$inclusion.p) == X2), 0.97)
  expect_equal(unname(fit$res_Gamma$beta), c(-0.8, 2), tolerance = 0.20)
})

test_that("listwise likelihood partition recovers clusters without imputing absent rows", {
  d <- simulate_early_paper(seed = 9206, N = 650)
  set.seed(77)
  missing_rows <- sample(seq_len(nrow(d$Z)), 0.2 * nrow(d$Z))
  Z_missing <- d$Z
  Z_missing[missing_rows, ] <- NA
  fit <- NULL
  invisible(capture.output(fit <- est_lucid(
    "early", d$G, Z_missing, d$Y_normal, CoY = d$CoY, K = 2,
    init_omic.data.model = "VVV", family = "normal", seed = 16,
    max_itr = 100, max_tot.itr = 300
  )))

  expect_gt(mean(max.col(fit$inclusion.p) == d$X), 0.94)
  expect_gt(mean(max.col(fit$inclusion.p[missing_rows, ]) == d$X[missing_rows]), 0.80)
  expect_true(all(is.na(fit$Z[missing_rows, ])))
  expect_equal(fit$res_Mu[, 1], c(-1.5, 1.5), tolerance = 0.25)
})

test_that("Equation 17 improves sporadic-value recovery under correlation", {
  set.seed(9207)
  N <- 500
  X <- sample(1:2, N, replace = TRUE)
  mu <- rbind(c(-1, -0.5), c(1, 0.5))
  covariance <- matrix(c(1, 0.8, 0.8, 1), 2)
  Z <- matrix(rnorm(N * 2), nrow = N) %*% chol(covariance) + mu[X, ]
  Z_missing <- Z
  Z_missing[, 2] <- NA
  initial <- Z_missing
  initial[, 2] <- mu[X, 2]
  p <- diag(2)[X, , drop = FALSE]
  sigma <- array(rep(covariance, 2), dim = c(2, 2, 2))
  imputed <- Istep_Z(initial, p, mu, sigma, !is.na(Z_missing), "parallel")
  rmse_equation17 <- sqrt(mean((imputed[, 2] - Z[, 2])^2))
  rmse_cluster_mean <- sqrt(mean((mu[X, 2] - Z[, 2])^2))

  expect_lt(rmse_equation17, 0.75 * rmse_cluster_mean)
})

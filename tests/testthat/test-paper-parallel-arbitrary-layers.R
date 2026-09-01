make_parallel_oracle_inputs <- function(K, N = 4) {
  G <- matrix(seq(-1, 1, length.out = N), ncol = 1)
  Beta <- lapply(seq_along(K), function(i) matrix(c(0.1 * i, -0.15), nrow = K[i] - 1))
  Mu <- lapply(seq_along(K), function(i) matrix(seq(-0.4, 0.4, length.out = K[i]), ncol = 1))
  Sigma <- lapply(K, function(k) array(rep(1, k), dim = c(1, 1, k)))
  Z <- lapply(seq_along(K), function(i) matrix(seq(-0.3, 0.5, length.out = N), ncol = 1))
  na_pattern <- lapply(Z, check_na)
  Delta <- parallel_delta_from_coef(c(-0.2, rep(0.1, sum(K - 1))), K, "gaussian", sd = 1)
  list(G = G, Beta = Beta, Mu = Mu, Sigma = Sigma, Z = Z,
       na_pattern = na_pattern, Delta = Delta, Y = matrix(seq(-1, 1, length.out = N)))
}

test_that("parallel E-step supports one layer without dropping model components", {
  K <- 2
  d <- make_parallel_oracle_inputs(K)
  got <- Estep(d$G, d$Z, d$Y, d$Beta, d$Mu, d$Sigma, d$Delta,
               "gaussian", TRUE, d$na_pattern)
  expected <- t(f_GtoX(d$G, d$Beta[[1]]) +
                  f_XtoZ(d$Z[[1]], d$Mu[[1]], d$Sigma[[1]])) +
    matrix(f_XtoY(d$Y, d$Delta, "gaussian"), nrow = 2)

  expect_equal(matrix(got, nrow = 2), expected)
  expect_equal(colSums(matrix(Estep_to_r(got, K, 4), 2, 4)), rep(1, 4))
})

test_that("parallel E-step and normalization support more than five layers", {
  K <- rep(2, 6)
  d <- make_parallel_oracle_inputs(K, N = 3)
  log_joint <- Estep(d$G, d$Z, d$Y, d$Beta, d$Mu, d$Sigma, d$Delta,
                     "gaussian", TRUE, d$na_pattern)
  r <- Estep_to_r(log_joint, K, 3)

  expect_equal(dim(log_joint), c(K, 3L))
  expect_true(all(is.finite(log_joint)))
  expect_equal(colSums(matrix(r, prod(K), 3)), rep(1, 3), tolerance = 1e-12)
  expect_equal(cal_loglik(log_joint),
               sum(vapply(1:3, function(i) LUCIDus:::safe_log_sum_exp(matrix(log_joint, prod(K), 3)[, i]),
                          numeric(1))))
})

test_that("g-computation joint probabilities have correct arbitrary-layer marginals", {
  K <- rep(2, 6)
  d <- make_parallel_oracle_inputs(K, N = 5)
  log_joint <- Estep_g(d$G, NULL, NULL, d$Beta, NULL, NULL, d$Delta,
                       "gaussian", FALSE, NULL)
  r <- Estep_to_r(log_joint, K, 5)

  for (layer in seq_along(K)) {
    expected <- exp(f_GtoX(d$G, d$Beta[[layer]]))
    expect_equal(compute_res_r(r, 5, layer), expected, tolerance = 1e-12)
  }
})

test_that("reference-coded state arrays support arbitrary layer counts", {
  K <- c(2, 3, 2, 2, 2, 2)
  coef <- c(0.7, seq_len(sum(K - 1)) / 10)
  got <- vec_to_array(K, coef)
  Delta <- parallel_delta_from_coef(coef, K, "gaussian")

  expect_equal(dim(got), K)
  expect_equal(got, Delta$eta)
})

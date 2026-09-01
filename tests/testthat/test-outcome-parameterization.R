# Test-only fixture helper: builds a flat non-reference-cluster design matrix
# from a joint responsibility array, purely to synthesize a plausible outcome
# below. This used to be a package-internal function (parallel_marginal_design())
# but had no call site anywhere in R/ -- only here -- so it was removed as
# dead code and the equivalent construction kept local to this test.
parallel_marginal_design <- function(r, K, N) {
  margins <- lapply(seq_along(K), function(layer) compute_res_r(r, N, layer))
  nonref <- lapply(seq_along(K), function(layer) {
    x <- margins[[layer]][, -1, drop = FALSE]
    colnames(x) <- paste0("Layer", layer, "_LC", seq.int(2, K[layer]))
    x
  })
  do.call(cbind, nonref)
}

test_that("early binary coefficients retain the fitted intercept and contrasts", {
  set.seed(101)
  n <- 150
  r <- matrix(rexp(n * 3), nrow = n)
  r <- r / rowSums(r)
  age <- rnorm(n)
  CoY <- matrix(age, ncol = 1, dimnames = list(NULL, "age"))
  Y <- rbinom(n, 1, plogis(-0.8 + 0.7 * r[, 2] - 0.3 * r[, 3] + 0.2 * age))

  observed <- fit_early_outcome(Y, r, CoY, 3, "age", "binary")
  long <- data.frame(
    Y = rep(Y, each = 3),
    state = factor(rep(1:3, times = n), levels = 1:3),
    age = rep(age, each = 3),
    w = as.vector(t(r))
  )
  expected <- suppressWarnings(glm(Y ~ state + age, family = binomial(),
                                    weights = w, data = long))

  expect_equal(unname(observed$beta), unname(coef(expected)), tolerance = 1e-10)
  expect_equal(observed$cluster_effect,
               unname(c(coef(expected)[1], coef(expected)[1] + coef(expected)[2:3])),
               tolerance = 1e-10)
  expect_equal(unname(observed$beta[2:3]), observed$cluster_effect[2:3] - observed$cluster_effect[1])
  expect_false(isTRUE(all.equal(exp(observed$beta[1]), exp(observed$beta[2]))))
})

test_that("early normal outcome uses absolute levels internally and contrasts in reports", {
  r <- matrix(c(.8, .2, .1, .9, .6, .4, .3, .7), ncol = 2, byrow = TRUE)
  Y <- c(-1, 0, 2, 3)
  gamma <- fit_early_outcome(Y, r, NULL, 2, NULL, "normal")
  expected <- colSums(r * Y) / colSums(r)

  expect_equal(gamma$cluster_effect, expected)
  expect_equal(unname(gamma$beta), c(expected[1], expected[2] - expected[1]))
  expect_equal(length(gamma$sigma), 2)
})

test_that("early relabeling preserves likelihood and all cluster-indexed parameters", {
  K <- 3
  N <- 7
  levels <- c(0.5, -1.2, -0.4)
  gamma <- early_gamma_from_levels(levels, covariate = 0.3,
                                   sigma = c(1.1, 0.8, 1.4), covariate_names = "age")
  beta <- rbind(c(0, 0), c(.2, -.1), c(-.4, .7))
  mu <- matrix(seq_len(9), nrow = K)
  sigma <- array(seq_len(27), dim = c(3, 3, K))
  CoY <- matrix(seq(-1, 1, length.out = N), ncol = 1)
  Y <- seq(-2, 2, length.out = N)

  before <- normal(K)$f.pYgX(Y, gamma, K, N, CoY, 1, itr = 2)
  out <- relabel_early_parameters(beta, mu, sigma, gamma, K)
  after <- normal(K)$f.pYgX(Y, out$gamma, K, N, CoY, 1, itr = 2)

  expect_identical(out$index, c(2L, 3L, 1L))
  expect_equal(after, before[, out$index])
  expect_equal(out$mu, mu[out$index, , drop = FALSE])
  expect_equal(unname(out$gamma$beta[-1]),
               c(out$gamma$cluster_effect[-1] - out$gamma$cluster_effect[1], 0.3))
  expect_equal(out$beta[1, ], c(0, 0))
})

test_that("parallel binary likelihood is Bernoulli by state and includes CoY", {
  K <- c(2, 2)
  Delta <- parallel_delta_from_coef(c(-1, 0.4, -0.3, 0.5), K,
                                    "binomial", covariate_names = "age")
  CoY <- matrix(c(-1, 1), ncol = 1)
  Y <- matrix(c(0, 1), ncol = 1)
  eta <- parallel_state_eta(Delta, CoY, N = 2)
  observed <- f_XtoY(Y, Delta, "binomial", CoY)
  expected <- array(dbinom(rep(c(0, 1), each = 4), 1,
                           plogis(as.numeric(eta)), log = TRUE), dim = c(K, 2))

  expect_equal(observed, expected)
  expect_gt(abs(sum(plogis(Delta$eta)) - 1), 1e-3)
  expect_equal(eta[, , 2] - eta[, , 1], matrix(1, 2, 2))
})

test_that("parallel normal likelihood uses additive state means and CoY", {
  K <- c(2, 3)
  Delta <- parallel_delta_from_coef(c(1, .4, -.2, .6, -.3), K,
                                    "gaussian", sd = .75,
                                    covariate_names = "age")
  CoY <- matrix(c(-2, .5, 1), ncol = 1)
  Y <- matrix(c(-1, 1, 2), ncol = 1)
  eta <- parallel_state_eta(Delta, CoY, N = 3)
  observed <- f_XtoY(Y, Delta, "gaussian", CoY)
  expected <- array(dnorm(rep(as.numeric(Y), each = prod(K)),
                          as.numeric(eta), .75, log = TRUE), dim = c(K, 3))

  expect_equal(observed, expected)
  expect_equal(eta[, , 3] - eta[, , 1], matrix(-.9, 2, 3))
})

test_that("parallel relabeling rebases effects without changing the model", {
  set.seed(102)
  K <- c(3, 2)
  N <- 5
  Delta <- parallel_delta_from_coef(c(-0.7, 0.9, -0.2, -0.5, 0.4), K,
                                    "gaussian", sd = 1.2,
                                    covariate_names = "age")
  Beta <- list(matrix(c(.2, -.1, .7, .3), nrow = 2), matrix(c(-.4, .6), nrow = 1))
  Mu <- list(matrix(seq_len(9), nrow = 3), matrix(seq_len(6), nrow = 2))
  Sigma <- list(array(seq_len(27), c(3, 3, 3)), array(seq_len(18), c(3, 3, 2)))
  r <- array(rexp(prod(K) * N), dim = c(K, N))
  for (n in seq_len(N)) r[, , n] <- r[, , n] / sum(r[, , n])

  out <- relabel_parallel_parameters(Beta, Mu, Sigma, Delta, r, K)
  expected_eta <- permute_parallel_array(Delta$eta, out$permutations)
  expected_r <- permute_parallel_array(r, out$permutations, include_observation = TRUE)

  expect_equal(out$Delta$eta, expected_eta)
  expect_equal(out$r, expected_r)
  expect_true(all(vapply(out$Delta$effects, function(x) all(diff(x) >= 0), logical(1))))
  for (i in seq_along(K)) {
    expect_equal(f_GtoX(matrix(c(-1, 0, 1), ncol = 1), out$Beta[[i]]),
                 f_GtoX(matrix(c(-1, 0, 1), ncol = 1), Beta[[i]])[, out$permutations[[i]]])
  }
})

test_that("parallel M-step maximizes the responsibility-weighted state likelihood", {
  set.seed(103)
  K <- c(2, 3)
  N <- 120
  r <- array(rexp(prod(K) * N), dim = c(K, N))
  for (n in seq_len(N)) r[, , n] <- r[, , n] / sum(r[, , n])
  CoY <- matrix(rnorm(N), ncol = 1, dimnames = list(NULL, "age"))
  design <- as.data.frame(parallel_marginal_design(r, K, N))
  Y <- rbinom(N, 1, plogis(-.4 + .5 * design[[1]] - .2 * design[[2]] +
                            .3 * design[[3]] + .2 * CoY[, 1]))

  observed <- Mstep_XtoY(Y, r, K, N, "binomial", CoY)
  grid <- expand.grid(1:2, 1:3)
  long <- data.frame(
    Y = rep(Y, each = 6),
    Layer1_LC2 = rep(as.numeric(grid[[1]] == 2), times = N),
    Layer2_LC2 = rep(as.numeric(grid[[2]] == 2), times = N),
    Layer2_LC3 = rep(as.numeric(grid[[2]] == 3), times = N),
    age = rep(CoY[, 1], each = 6),
    w = as.vector(matrix(r, 6, N))
  )
  expected <- suppressWarnings(glm(Y ~ Layer1_LC2 + Layer2_LC2 + Layer2_LC3 + age,
                                    data = long, weights = w, family = binomial()))
  expect_equal(unname(parallel_delta_coef(observed$Gamma)),
               unname(coef(expected)), tolerance = 1e-10)
})

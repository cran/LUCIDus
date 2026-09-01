test_that("early binary M-step equals a weighted latent-state GLM", {
  set.seed(8101)
  N <- 180
  K <- 3
  r <- matrix(rexp(N * K), N, K)
  r <- r / rowSums(r)
  age <- rnorm(N)
  CoY <- matrix(age, ncol = 1, dimnames = list(NULL, "age"))
  Y <- rbinom(N, 1, plogis(-0.7 + 0.4 * age))

  got <- fit_early_outcome(Y, r, CoY, K, "age", "binary")
  long <- data.frame(
    Y = rep(Y, each = K),
    state = factor(rep(seq_len(K), times = N), levels = seq_len(K)),
    age = rep(age, each = K),
    w = as.vector(t(r))
  )
  oracle <- suppressWarnings(glm(Y ~ state + age, data = long, weights = w,
                                 family = binomial()))
  posterior_score <- glm(Y ~ r[, 2] + r[, 3] + age, family = binomial())

  expect_equal(unname(got$beta), unname(coef(oracle)), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(unname(got$beta), unname(coef(posterior_score)))))
})

test_that("early Gaussian M-step recovers cluster means, CoY effect, and variances", {
  set.seed(8102)
  N <- 2400
  state <- sample(1:3, N, replace = TRUE, prob = c(.3, .45, .25))
  r <- diag(3)[state, , drop = FALSE]
  age <- rnorm(N)
  levels <- c(-1.2, 0.4, 1.7)
  alpha <- 0.65
  sigma <- c(0.5, 1.1, 1.8)
  Y <- levels[state] + alpha * age + rnorm(N, sd = sigma[state])
  CoY <- matrix(age, ncol = 1, dimnames = list(NULL, "age"))

  got <- fit_early_outcome(Y, r, CoY, 3, "age", "normal")

  expect_equal(got$cluster_effect, levels, tolerance = 0.12)
  expect_equal(got$covariate, alpha, tolerance = 0.06)
  expect_equal(unname(got$sigma), sigma, tolerance = 0.10)
  expect_equal(unname(got$beta[2:3]), got$cluster_effect[2:3] - got$cluster_effect[1],
               tolerance = 1e-12)
})

test_that("parallel Gaussian M-step equals weighted additive state regression", {
  set.seed(8103)
  K <- c(2, 3)
  N <- 160
  S <- prod(K)
  r <- array(rexp(S * N), dim = c(K, N))
  for (i in seq_len(N)) r[, , i] <- r[, , i] / sum(r[, , i])
  age <- rnorm(N)
  CoY <- matrix(age, ncol = 1, dimnames = list(NULL, "age"))
  Y <- rnorm(N, 0.3 + 0.5 * age, 1.2)

  got <- fit_parallel_outcome(Y, r, K, N, "gaussian", CoY)$Gamma
  grid <- expand.grid(1:2, 1:3)
  long <- data.frame(
    Y = rep(Y, each = S),
    Layer1_LC2 = rep(as.numeric(grid[[1]] == 2), times = N),
    Layer2_LC2 = rep(as.numeric(grid[[2]] == 2), times = N),
    Layer2_LC3 = rep(as.numeric(grid[[2]] == 3), times = N),
    age = rep(age, each = S),
    w = as.vector(matrix(r, S, N))
  )
  oracle <- lm(Y ~ Layer1_LC2 + Layer2_LC2 + Layer2_LC3 + age,
               data = long, weights = w)
  oracle_sd <- sqrt(sum(long$w * residuals(oracle)^2) / sum(long$w))

  expect_equal(unname(got$coef), unname(coef(oracle)), tolerance = 1e-10)
  expect_equal(got$sd, oracle_sd, tolerance = 1e-12)
})

test_that("CoY contributes to the early likelihood from the first E-step", {
  K <- 2
  gamma <- early_gamma_from_levels(c(-0.4, 0.8), covariate = 0.7,
                                   sigma = c(0.6, 1.2), covariate_names = "age")
  Y <- c(-1, 0.5, 2)
  CoY <- matrix(c(-1, 0, 2), ncol = 1)
  got <- normal(K)$f.pYgX(Y, gamma, K, 3, CoY, 1, itr = 1)
  expected <- vapply(1:K, function(k) {
    dnorm(Y, gamma$cluster_effect[k] + 0.7 * CoY[, 1], gamma$sigma[k], log = TRUE)
  }, numeric(3))

  expect_equal(got, expected)
})

test_that("early EM does not repeatedly restart after max_itr without numerical failure", {
  set.seed(7301)
  N <- 30
  G <- matrix(rnorm(N * 2), nrow = N)
  Z <- matrix(rnorm(N * 4), nrow = N)
  Y <- rnorm(N)

  out <- suppressWarnings(capture.output(
    fit <- estimate_lucid(
      lucid_model = "early",
      G = G,
      Z = Z,
      Y = Y,
      family = "normal",
      K = 2,
      max_itr = 1,
      max_tot.itr = 50,
      tol = 1e-12,
      seed = 7301,
      verbose = TRUE
    )
  ))

  n_restart_msgs <- sum(grepl("^Fitting Early Integration LUCID model", out))
  expect_equal(n_restart_msgs, 1)
  expect_s3_class(fit, "early_lucid")
})

test_that("parallel EM does not repeatedly restart after max_itr without numerical failure", {
  set.seed(7302)
  N <- 30
  G <- matrix(rnorm(N * 2), nrow = N)
  Z1 <- matrix(rnorm(N * 3), nrow = N)
  Z2 <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  out <- suppressWarnings(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      family = "normal",
      K = c(2, 2),
      max_itr = 1,
      max_tot.itr = 50,
      tol = 1e-12,
      seed = 7302,
      verbose = TRUE
    )
  ))

  n_restart_msgs <- sum(grepl("^Fitting LUCID in Parallel model", out))
  expect_equal(n_restart_msgs, 1)
  expect_s3_class(fit, "lucid_parallel")
})

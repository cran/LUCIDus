test_that("early verbose FALSE prints concise start/finish and penalty selection counts", {
  set.seed(8601)
  N <- 40
  G <- matrix(rnorm(N * 3), nrow = N)
  Z <- matrix(rnorm(N * 5), nrow = N)
  Y <- rnorm(N)

  out <- suppressWarnings(capture.output(
    fit <- estimate_lucid(
      lucid_model = "early",
      G = G,
      Z = Z,
      Y = Y,
      family = "normal",
      K = 2,
      Rho_G = 0.01,
      max_itr = 6,
      max_tot.itr = 20,
      tol = 1e-2,
      seed = 8601,
      verbose = FALSE
    )
  ))

  expect_s3_class(fit, "early_lucid")
  expect_true(any(grepl("^Fitting LUCID early model", out)))
  expect_true(any(grepl("^Finished LUCID early model", out)))
  expect_true(any(grepl("Selected G:", out)))
  expect_false(any(grepl("^iteration\\s+[0-9]+", out)))
  expect_false(any(grepl("^\\.+$", out)))
})

test_that("parallel verbose TRUE prints iteration log-likelihood updates", {
  set.seed(8602)
  N <- 36
  G <- matrix(rnorm(N * 3), nrow = N)
  Z <- list(matrix(rnorm(N * 4), nrow = N), matrix(rnorm(N * 4), nrow = N))
  Y <- rnorm(N)

  out <- suppressWarnings(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = Z,
      Y = Y,
      family = "normal",
      K = c(2, 2),
      max_itr = 4,
      max_tot.itr = 20,
      tol = 1e-2,
      seed = 8602,
      verbose = TRUE
    )
  ))

  expect_s3_class(fit, "lucid_parallel")
  expect_true(any(grepl("^Fitting LUCID in Parallel model", out)))
  expect_true(any(grepl("^iteration\\s+[0-9]+: log-likelihood", out)))
})

test_that("early verbose TRUE prints iteration log-likelihood updates, matching parallel", {
  # Early used to print "iteration N: E-step finished." and only report the
  # log-likelihood inside a differently-worded M-step line ("loglike"/
  # "penalized loglike"). This pins the fix: the per-iteration line now uses
  # the same "log-likelihood" wording parallel's does.
  set.seed(8604)
  N <- 36
  G <- matrix(rnorm(N * 3), nrow = N)
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
      max_itr = 4,
      max_tot.itr = 20,
      tol = 1e-2,
      seed = 8604,
      verbose = TRUE
    )
  ))

  expect_s3_class(fit, "early_lucid")
  expect_true(any(grepl("^iteration\\s+[0-9]+: log-likelihood", out)))
})

test_that("serial verbose TRUE uses Stage naming, not Sub Model, in the pre-stage banner", {
  set.seed(8605)
  N <- 36
  G <- matrix(rnorm(N * 3), nrow = N)
  Z <- list(matrix(rnorm(N * 4), nrow = N), matrix(rnorm(N * 4), nrow = N))
  Y <- rnorm(N)

  out <- suppressWarnings(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = Z,
      Y = Y,
      family = "normal",
      K = list(2, 2),
      max_itr = 4,
      max_tot.itr = 20,
      tol = 1e-2,
      seed = 8605,
      verbose = TRUE
    )
  ))

  expect_s3_class(fit, "lucid_serial")
  expect_true(any(grepl("^Fitting LUCID serial model \\(Stage 1/2\\)", out)))
  expect_true(any(grepl("^Fitting LUCID serial model \\(Stage 2/2\\)", out)))
  expect_false(any(grepl("Sub Model", out)))
})

test_that("serial verbose FALSE prints concise stage-level progress", {
  set.seed(8603)
  N <- 36
  G <- matrix(rnorm(N * 3), nrow = N)
  Z <- list(
    matrix(rnorm(N * 4), nrow = N),
    matrix(rnorm(N * 4), nrow = N)
  )
  Y <- rnorm(N)

  out <- suppressWarnings(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = Z,
      Y = Y,
      family = "normal",
      K = list(2, 2),
      max_itr = 5,
      max_tot.itr = 20,
      tol = 1e-2,
      seed = 8603,
      verbose = FALSE
    )
  ))

  expect_s3_class(fit, "lucid_serial")
  expect_true(any(grepl("^Fitting LUCID serial model", out)))
  expect_true(any(grepl("^  Stage 1/2 \\(early\\) finished: log-likelihood = ", out)))
  expect_true(any(grepl("^  Stage 2/2 \\(early\\) finished: log-likelihood = ", out)))
  expect_true(any(grepl("^Finished LUCID serial model\\.", out)))
  expect_false(any(grepl("^iteration\\s+[0-9]+", out)))
})

# Focused tests for summary missing-data reporting

test_that("early summary stores and prints missing profile", {
  set.seed(1008)
  G <- matrix(rnorm(120), nrow = 30)
  Z <- matrix(rnorm(240), nrow = 30)
  Y <- rnorm(30)
  Z[1, ] <- NA
  Z[2:3, 1] <- NA

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008,
      max_itr = 10,
      tol = 1e-1
    )
  )))

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_early")
  expect_equal(s$missing_data$listwise_rows, 1)
  expect_equal(s$missing_data$sporadic_rows, 2)
  txt <- capture.output(print(s))
  expect_true(any(grepl("Missing-data profile", txt)))
})

test_that("parallel summary stores and prints per-layer missing profile", {
  set.seed(1008)
  G <- matrix(rnorm(120), nrow = 30)
  Z1 <- matrix(rnorm(180), nrow = 30)
  Z2 <- matrix(rnorm(180), nrow = 30)
  Y <- rnorm(30)
  Z1[1, ] <- NA
  Z1[2, 1] <- NA
  Z2[3, ] <- NA
  Z2[4:5, 1] <- NA

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = c(2, 2),
      family = "normal",
      init_impute = "mix",
      seed = 1008,
      max_itr = 8,
      tol = 1e-1
    )
  )))

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_parallel")
  expect_equal(s$missing_data$layer_summary$listwise_rows[1], 1)
  expect_equal(s$missing_data$layer_summary$sporadic_rows[1], 1)
  expect_equal(s$missing_data$layer_summary$listwise_rows[2], 1)
  expect_equal(s$missing_data$layer_summary$sporadic_rows[2], 2)
  txt <- capture.output(print(s))
  expect_true(any(grepl("Missing-data profile by layer", txt)))
})

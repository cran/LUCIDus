# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# summary()'s printed bootstrap CI tables (G-to-X, cluster-to-Y, and
# cluster-specific omics means) mark a row "*" when its normal-theory CI
# excludes 0. The column only exists when boot.se is supplied -- there is no
# CI to test otherwise. sig is always determined from the normal-theory
# interval alone, so the percentile-interval columns are dropped from the
# printed table entirely -- showing them next to a sig column that never
# used them would suggest they mattered when they don't.

test_that(".append_significance adds '*' only where the interval excludes 0", {
  df <- data.frame(estimate = c(1, 0.1, -1),
                   norm_lower = c(0.5, -0.2, -1.5),
                   norm_upper = c(1.5, 0.4, -0.5))
  out <- LUCIDus:::.append_significance(df)
  expect_equal(out$sig, c("*", "", "*"))
})

test_that(".append_significance drops the percentile-interval columns", {
  df <- data.frame(estimate = c(1, -1),
                   norm_lower = c(0.5, -1.5), norm_upper = c(1.5, -0.5),
                   perc_lower = c(0.4, -1.6), perc_upper = c(1.6, -0.4))
  out <- LUCIDus:::.append_significance(df)
  expect_false(any(c("perc_lower", "perc_upper") %in% colnames(out)))
  expect_equal(out$sig, c("*", "*"))
  expect_equal(colnames(out), c("estimate", "norm_lower", "norm_upper", "sig"))
})

test_that(".append_significance is a no-op without CI columns, and is matrix-safe", {
  df <- data.frame(estimate = c(1, 2))
  expect_identical(LUCIDus:::.append_significance(df), df)

  m <- matrix(c(1, 2, 0.5, -0.5, 1.5, 0.5), nrow = 2,
             dimnames = list(c("a", "b"), c("estimate", "norm_lower", "norm_upper")))
  out <- LUCIDus:::.append_significance(m)
  expect_s3_class(out, "data.frame")
  expect_equal(out$sig, c("*", ""))
})

test_that(".append_significance treats a non-finite bound as not significant", {
  df <- data.frame(estimate = 1, norm_lower = NA_real_, norm_upper = 2)
  out <- LUCIDus:::.append_significance(df)
  expect_equal(out$sig, "")
})

test_that("summary() shows no sig column without boot.se, and a sig column with it", {
  set.seed(71)
  N <- 200
  G <- matrix(rnorm(N * 2), nrow = N); colnames(G) <- c("g1", "g2")
  X <- rbinom(N, 1, plogis(2 * G[, 1])) + 1
  Z <- matrix(rnorm(N * 2, c(-2, 2)[X]), nrow = N); colnames(Z) <- c("z1", "z2")
  Y <- rnorm(N, c(-1, 1)[X])
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(G = G, Z = Z, Y = Y, K = 2, lucid_model = "early",
                          family = "normal", init_omic.data.model = "VVV", seed = 71)
  ))))

  out_no_ci <- capture.output(summary(fit))
  expect_false(any(grepl("\\bsig\\b", out_no_ci)))

  boot <- suppressWarnings(boot_lucid(G = G, Z = Z, Y = Y, model = fit, R = 30))
  out_ci <- capture.output(summary(fit, boot.se = boot))
  # the well-separated data should produce at least one significant row in
  # each of the three tables the sig column now appears in
  sig_lines <- grep("\\*\\s*$", out_ci, value = TRUE)
  expect_true(length(sig_lines) > 0)
  expect_true(any(grepl("\\bsig\\b", out_ci)))
  expect_false(any(grepl("perc_lower|perc_upper", out_ci)))
})

test_that("summary() shows a sig column for the parallel model's per-layer tables", {
  set.seed(72)
  N <- 200
  G <- matrix(rnorm(N * 2), nrow = N)
  Z1 <- matrix(rnorm(N * 2), nrow = N)
  Z2 <- matrix(rnorm(N * 2), nrow = N)
  Y <- rnorm(N)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(G = G, Z = list(Z1, Z2), Y = Y, K = c(2, 2),
                          lucid_model = "parallel", family = "normal",
                          init_omic.data.model = "VVV", seed = 72)
  ))))
  boot <- suppressWarnings(boot_lucid(G = G, Z = list(Z1, Z2), Y = Y, model = fit, R = 15))
  out <- capture.output(summary(fit, boot.se = boot))
  expect_true(any(grepl("\\bsig\\b", out)))
})

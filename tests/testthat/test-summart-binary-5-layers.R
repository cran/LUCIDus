# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# LUCID summary smoke test for early model (binary)

test_that("summary works for tuned early binary model", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z_e <- matrix(rnorm(1000), nrow = 100)
  Y <- rbinom(n = 100, size = 1, prob = 0.25)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)

  suppressWarnings(invisible(capture.output(
    fit1 <- lucid(
      G = G, Z = Z_e, Y = Y, K = 2:3,
      CoG = CoG, CoY = CoY,
      lucid_model = "early",
      family = "binary",
      init_omic.data.model = "VVV",
      seed = i,
      useY = TRUE
    )
  )))

  sum_fit1 <- summary(fit1)
  expect_s3_class(fit1, "early_lucid")
  expect_s3_class(sum_fit1, "sumlucid_early")
  expect_true(is.finite(sum_fit1$model_fit$BIC))
  expect_equal(nrow(fit1$inclusion.p), nrow(G))
})

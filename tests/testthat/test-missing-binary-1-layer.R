# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# LUCID - 1 omics, binary outcome with missing data

test_that("early binary handles missing Z with stable posterior output", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z <- matrix(rnorm(1000), nrow = 100)
  miss_idx <- sample(seq_len(length(Z)), 40, replace = FALSE)
  Z[miss_idx] <- NA
  Z[2, ] <- NA
  Y <- rbinom(n = 100, size = 1, prob = 0.25)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)

  suppressWarnings(invisible(capture.output(
    fit1 <- estimate_lucid(
      G = G, Z = Z, Y = Y, K = 2,
      CoG = CoG, CoY = CoY,
      lucid_model = "early",
      family = "binary",
      seed = i,
      useY = TRUE
    )
  )))

  expect_s3_class(fit1, "early_lucid")
  expect_equal(dim(fit1$inclusion.p), c(nrow(G), 2))
  expect_equal(rowSums(fit1$inclusion.p), rep(1, nrow(G)), tolerance = 1e-6)
  expect_true(all(is.finite(fit1$inclusion.p)))
  expect_equal(sum(is.na(fit1$Z)), ncol(Z))
  expect_true(all(is.na(fit1$Z[2, ])))
  expect_true(all(is.finite(fit1$Z[-2, ])))
})

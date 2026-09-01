# Heavy: fits multiple LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# K >= 2 was already enforced when fitting through lucid() (early), through
# tune_lucid()'s parallel branch, and through estimate_lucid()'s serial
# branch -- but est_lucid() itself (early and parallel) and tune_lucid_auxi()'s
# own early branch never checked it, so calling estimate_lucid()/tune_lucid()
# directly with K = 1 (which the docs explicitly support) reached deep EM
# machinery before failing with an obscure error instead of a clear one.

test_that("estimate_lucid() rejects K = 1 for early with a clear message", {
  set.seed(101)
  N <- 60
  G <- matrix(rnorm(N * 2), nrow = N)
  Z <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  expect_error(
    estimate_lucid(lucid_model = "early", G = G, Z = Z, Y = Y,
                   family = "normal", K = 1, init_omic.data.model = "VVV"),
    "greater or equal than 2"
  )
})

test_that("estimate_lucid() rejects K = 1 for parallel with a clear message", {
  set.seed(102)
  N <- 60
  G <- matrix(rnorm(N * 2), nrow = N)
  Z <- list(matrix(rnorm(N * 3), nrow = N), matrix(rnorm(N * 3), nrow = N))
  Y <- rnorm(N)

  expect_error(
    estimate_lucid(lucid_model = "parallel", G = G, Z = Z, Y = Y,
                   family = "normal", K = c(1, 2), init_omic.data.model = "VVV"),
    "greater or equal than 2"
  )
})

test_that("tune_lucid() rejects a K grid containing 1 for early, up front", {
  set.seed(103)
  N <- 60
  G <- matrix(rnorm(N * 2), nrow = N)
  Z <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  expect_error(
    suppressWarnings(suppressMessages(invisible(capture.output(
      tune_lucid(G = G, Z = Z, Y = Y, lucid_model = "early", K = 1:3)
    )))),
    "greater or equal than 2"
  )
})

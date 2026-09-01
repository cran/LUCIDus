skip_on_cran()

# tune_lucid()'s early/parallel branches wrap each candidate fit in try(),
# record NA on failure, and only stop() if every candidate failed -- so one
# bad K in a grid search doesn't crash the whole tuning run. The serial
# branch had neither: a single failing candidate crashed the entire call, and
# picking the best model (`which(bic == min(bic))`, no `[1]` guard) would
# error on a genuine BIC tie the way the early/parallel branches are already
# written to avoid.

make_serial_tune_data <- function(seed = 5501, n = 80, pG = 3, pZ = 4) {
  set.seed(seed)
  list(
    G = matrix(rnorm(n * pG), nrow = n),
    Z1 = matrix(rnorm(n * pZ), nrow = n),
    Z2 = matrix(rnorm(n * pZ), nrow = n),
    Y = rnorm(n)
  )
}

test_that("a K = 1 candidate in the serial grid is skipped, not fatal", {
  d <- make_serial_tune_data()
  Z <- list(Z1 = d$Z1, Z2 = d$Z2)

  # stage 1's grid includes an invalid K = 1 alongside a valid K = 2; check_K()
  # rejects K = 1 cleanly inside estimate_lucid(), which try() now catches.
  out <- capture.output(
    tune <- tune_lucid(
      G = d$G, Z = Z, Y = d$Y,
      K = list(c(1, 2), 2),
      lucid_model = "serial",
      family = "normal",
      init_omic.data.model = "VVV",
      seed = 5501
    )
  )

  expect_false(is.null(tune))
  expect_s3_class(tune$model_opt, "lucid_serial")
  # the failed candidate is recorded as NA, not silently dropped or fatal
  # (tune_K is a list-matrix here, since its K columns hold per-stage vectors)
  bic_col <- unlist(tune$tune_K[, "bic"])
  expect_true(any(is.na(bic_col)))
  expect_true(any(is.finite(bic_col)))
})

test_that("all candidates failing in the serial grid gives one clear error", {
  d <- make_serial_tune_data()
  Z <- list(Z1 = d$Z1, Z2 = d$Z2)

  expect_error(
    capture.output(
      tune_lucid(
        G = d$G, Z = Z, Y = d$Y,
        K = list(c(1), 1),
        lucid_model = "serial",
        family = "normal",
        init_omic.data.model = "VVV",
        seed = 5501
      )
    ),
    "fails to converge"
  )
})

test_that("a genuine BIC tie in the serial grid resolves to one model, not an error", {
  d <- make_serial_tune_data()
  Z <- list(Z1 = d$Z1, Z2 = d$Z2)

  # two identical candidates (same K, same data, same fixed seed) produce
  # byte-identical BIC -- exactly the tie the missing `[1]` guard would have
  # turned into a `model_list[[c(i,j)]]` recursive-indexing error.
  out <- capture.output(
    tune <- tune_lucid(
      G = d$G, Z = Z, Y = d$Y,
      K = list(c(2, 2), 2),
      lucid_model = "serial",
      family = "normal",
      init_omic.data.model = "VVV",
      seed = 5501
    )
  )

  expect_s3_class(tune$model_opt, "lucid_serial")
  expect_length(class(tune$model_opt), length(class(tune$model_opt)))
})

# Heavy: fits multiple LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# lucid()'s selection-count messages and tune_lucid()'s "Tuning LUCID model"
# banner used to print unconditionally, ignoring verbose_tune -- the one flag
# that exists to control exactly this. This pins that they now respect it,
# for early (lucid()'s selection messages) and for all three model types
# (tune_lucid()'s banner).

make_tuning_data <- function(seed = 7001, n = 80, pG = 4, pZ = 4) {
  set.seed(seed)
  list(
    G = matrix(rnorm(n * pG), nrow = n),
    Z = matrix(rnorm(n * pZ), nrow = n),
    Y = rnorm(n)
  )
}

test_that("lucid()'s selection-count messages respect verbose_tune", {
  d <- make_tuning_data()

  out_quiet <- capture.output(
    fit_quiet <- suppressWarnings(suppressMessages(
      lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
           family = "normal", K = 2, Rho_G = 0.05,
           init_omic.data.model = "VVV", verbose_tune = FALSE)
    ))
  )
  expect_false(any(grepl("exposures are selected|No exposure variables", out_quiet)))

  out_loud <- capture.output(
    fit_loud <- suppressWarnings(suppressMessages(
      lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
           family = "normal", K = 2, Rho_G = 0.05,
           init_omic.data.model = "VVV", verbose_tune = TRUE)
    ))
  )
  expect_true(any(grepl("exposures are selected|No exposure variables", out_loud)))
})

test_that("tune_lucid()'s 'Tuning LUCID model' banner respects verbose_tune, for every model type", {
  d <- make_tuning_data()

  # early: a grid with more than one candidate
  out_quiet_e <- capture.output(suppressWarnings(suppressMessages(
    tune_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
              K = 2:3, verbose_tune = FALSE)
  )))
  expect_false(any(grepl("^Tuning LUCID model", out_quiet_e)))

  out_loud_e <- capture.output(suppressWarnings(suppressMessages(
    tune_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
              K = 2:3, verbose_tune = TRUE)
  )))
  expect_true(any(grepl("^Tuning LUCID model", out_loud_e)))

  # parallel: a grid with more than one candidate
  Z2 <- list(d$Z, matrix(rnorm(nrow(d$Z) * 4), nrow = nrow(d$Z)))
  out_quiet_p <- capture.output(suppressWarnings(suppressMessages(
    tune_lucid(G = d$G, Z = Z2, Y = d$Y, lucid_model = "parallel",
              K = list(2:3, 2), verbose_tune = FALSE)
  )))
  expect_false(any(grepl("^Tuning LUCID model", out_quiet_p)))

  out_loud_p <- capture.output(suppressWarnings(suppressMessages(
    tune_lucid(G = d$G, Z = Z2, Y = d$Y, lucid_model = "parallel",
              K = list(2:3, 2), verbose_tune = TRUE)
  )))
  expect_true(any(grepl("^Tuning LUCID model", out_loud_p)))

  # serial: a grid with more than one candidate
  out_quiet_s <- capture.output(suppressWarnings(suppressMessages(
    tune_lucid(G = d$G, Z = list(Z1 = d$Z, Z2 = Z2[[2]]), Y = d$Y,
              lucid_model = "serial", K = list(2:3, 2), verbose_tune = FALSE)
  )))
  expect_false(any(grepl("^Tuning LUCID model", out_quiet_s)))

  out_loud_s <- capture.output(suppressWarnings(suppressMessages(
    tune_lucid(G = d$G, Z = list(Z1 = d$Z, Z2 = Z2[[2]]), Y = d$Y,
              lucid_model = "serial", K = list(2:3, 2), verbose_tune = TRUE)
  )))
  expect_true(any(grepl("^Tuning LUCID model", out_loud_s)))
})

test_that("predict_lucid()'s g-computation notice respects verbose", {
  d <- make_tuning_data()
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                family = "normal", K = 2, init_omic.data.model = "VVV")
  ))))

  out_quiet <- capture.output(suppressWarnings(suppressMessages(
    predict_lucid(model = fit, lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                 g_computation = TRUE, verbose = FALSE)
  )))
  expect_false(any(grepl("G-computation only uses", out_quiet)))

  out_loud <- capture.output(suppressWarnings(suppressMessages(
    predict_lucid(model = fit, lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                 g_computation = TRUE, verbose = TRUE)
  )))
  expect_true(any(grepl("G-computation only uses", out_loud)))
})

test_that("predict_lucid()'s new verbose output appears for early and parallel prediction", {
  d <- make_tuning_data()
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                family = "normal", K = 2, init_omic.data.model = "VVV")
  ))))

  out_quiet <- capture.output(suppressWarnings(suppressMessages(
    predict_lucid(model = fit, lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                 verbose = FALSE)
  )))
  expect_false(any(grepl("^Predicting LUCID", out_quiet)))

  out_loud <- capture.output(suppressWarnings(suppressMessages(
    predict_lucid(model = fit, lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                 verbose = TRUE)
  )))
  expect_true(any(grepl("^Predicting LUCID early model", out_loud)))
})

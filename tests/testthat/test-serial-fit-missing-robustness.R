# Focused, fast regression tests for serial fitting + missing-data flow.

make_serial_smoke_data <- function(n = 40, pG = 5, pZ = 4, seed = 2026) {
  set.seed(seed)
  G <- matrix(rnorm(n * pG), nrow = n)
  Z1 <- matrix(rnorm(n * pZ), nrow = n)
  Z2 <- matrix(rnorm(n * pZ), nrow = n)
  Z3 <- matrix(rnorm(n * pZ), nrow = n)
  Z4 <- matrix(rnorm(n * pZ), nrow = n)
  Z1[1, ] <- NA
  Z2[2, 1] <- NA
  Z3[3, ] <- NA
  Yn <- rnorm(n)
  Yb <- rbinom(n, 1, 0.4)
  CoG <- matrix(rnorm(n * 2), nrow = n)
  CoY <- matrix(rnorm(n * 2), nrow = n)
  list(G = G, Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Yn = Yn, Yb = Yb, CoG = CoG, CoY = CoY)
}

test_that("serial estimate_lucid fits mixed topology with missing data and stores stage summaries", {
  d <- make_serial_smoke_data(seed = 1008)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = d$G,
      Z = list(d$Z1, list(d$Z2, d$Z3), d$Z4),
      Y = d$Yb,
      CoG = d$CoG,
      CoY = d$CoY,
      K = list(2, list(2, 2), 2),
      family = "binary",
      max_itr = 5,
      max_tot.itr = 20,
      seed = 1008
    )
  )))

  expect_s3_class(fit, "lucid_serial")
  expect_equal(length(fit$submodel), 3)
  expect_true(is.list(fit$missing_summary))
  expect_equal(fit$missing_summary$n_stages, 3)
  expect_equal(length(fit$missing_summary$stage), 3)
  expect_true(all(sapply(fit$submodel, function(m) !is.null(m$missing_summary))))
})

test_that("serial accepts data.frame omics blocks when K topology is valid", {
  d <- make_serial_smoke_data(seed = 1010)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = d$G,
      Z = list(as.data.frame(d$Z1), list(as.data.frame(d$Z2), as.data.frame(d$Z3))),
      Y = d$Yn,
      CoG = d$CoG,
      CoY = d$CoY,
      K = list(2, list(2, 2)),
      family = "normal",
      max_itr = 5,
      max_tot.itr = 20,
      seed = 1010
    )
  )))

  expect_s3_class(fit, "lucid_serial")
  expect_equal(length(fit$submodel), 2)
})

test_that("serial wrapper tune_lucid/lucid works and handles stage with single-column G under Rho_G", {
  d <- make_serial_smoke_data(seed = 1012)

  suppressWarnings(invisible(capture.output(
    tune <- tune_lucid(
      G = d$G,
      Z = list(d$Z1, d$Z2),
      Y = d$Yn,
      CoG = d$CoG,
      CoY = d$CoY,
      family = "normal",
      K = list(2, 2:3),
      lucid_model = "serial",
      Rho_G = 0.2,
      max_itr = 4,
      max_tot.itr = 16,
      seed = 1012
    )
  )))

  suppressWarnings(invisible(capture.output(
    fit <- lucid(
      G = d$G,
      Z = list(d$Z1, d$Z2),
      Y = d$Yn,
      CoG = d$CoG,
      CoY = d$CoY,
      family = "normal",
      K = list(2, 2:3),
      lucid_model = "serial",
      Rho_G = 0.2,
      max_itr = 4,
      max_tot.itr = 16,
      seed = 1012
    )
  )))

  expect_true(nrow(as.data.frame(tune$tune_K)) >= 1)
  expect_s3_class(fit, "lucid_serial")
  expect_equal(fit$submodel[[1]]$Rho$Rho_G, 0.2)
  expect_equal(fit$submodel[[2]]$Rho$Rho_G, 0)
})

test_that("serial summary uses structured style and keeps detailed stage parameter tables", {
  d <- make_serial_smoke_data(seed = 1014)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = d$G,
      Z = list(d$Z1, list(d$Z2, d$Z3), d$Z4),
      Y = d$Yn,
      CoG = d$CoG,
      CoY = d$CoY,
      K = list(2, list(2, 2), 2),
      family = "normal",
      max_itr = 4,
      max_tot.itr = 16,
      seed = 1014
    )
  )))

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_serial")
  expect_true(is.list(s$model_info))
  expect_equal(s$model_info$n_stages, 3)
  expect_true(is.list(s$stage_summary))
  expect_equal(length(s$stage_summary), 3)
  expect_true(any(sapply(s$stage_summary, inherits, what = "sumlucid_parallel")))
  expect_true(any(sapply(s$stage_summary, inherits, what = "sumlucid_early")))

  txt <- capture.output(print(s))
  expect_true(any(grepl("LUCID Serial: Model Summary", txt)))
  expect_true(any(grepl("Missing-data profile by stage", txt)))
  expect_true(any(grepl("Stage-wise detailed parameter estimates", txt)))
  expect_true(any(grepl("LUCID Early Integration: Model Summary", txt)))
  expect_true(any(grepl("LUCID Parallel: Model Summary", txt)))
})

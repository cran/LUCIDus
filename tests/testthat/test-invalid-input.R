# Crazy values in the inputs must be rejected with a message that says what is
# wrong, not absorbed silently.
#
# LUCID models missingness in the OMICS data only. Missing or non-finite values
# anywhere else -- exposures, outcome, covariates -- are not handled by any code
# path, and the symptom used to be either a fit that ran for minutes without
# converging or one that returned confidently wrong estimates with nothing in
# the output to indicate a problem. A single NA in Y was enough to do it.

make_data <- function(model = "early", n_layer = 1L, N = 300L, ...) {
  sim_lucid(model, N = N, K = 2, M = 4, n_layer = n_layer, seed = 900, ...)
}

fit_it <- function(d, model = "early", K = 2, family = "normal", ...) {
  suppressWarnings(suppressMessages(invisible(capture.output(
    out <- estimate_lucid(lucid_model = model, G = d$G, Z = d$Z, Y = d$Y,
                          family = family, K = K,
                          init_omic.data.model = "VVV", ...)
  ))))
  out
}

test_that("a missing value in the outcome is rejected", {
  d <- make_data()
  d$Y[3] <- NA
  expect_error(fit_it(d), "Y contains 1 missing value")
})

test_that("a missing value in the exposures is rejected", {
  d <- make_data()
  d$G[5, 1] <- NA
  expect_error(fit_it(d), "G contains 1 missing value")
})

test_that("missing values in covariates are rejected on either arm", {
  d <- make_data(n_CoG = 2, n_CoY = 2)
  bad_cog <- d$CoG; bad_cog[2, 1] <- NA
  expect_error(fit_it(d, CoG = bad_cog, CoY = d$CoY), "CoG contains 1 missing value")

  bad_coy <- d$CoY; bad_coy[4, 2] <- NA
  expect_error(fit_it(d, CoG = d$CoG, CoY = bad_coy), "CoY contains 1 missing value")
})

test_that("non-finite values are rejected wherever they appear", {
  d <- make_data()
  d$G[7, 1] <- Inf
  expect_error(fit_it(d), "G contains non-finite values")

  d2 <- make_data()
  d2$Y[9] <- -Inf
  expect_error(fit_it(d2), "Y contains non-finite values")

  d3 <- make_data()
  d3$G[2, 2] <- NaN
  expect_error(fit_it(d3), "G contains 1 missing value")   # NaN counts as NA
})

test_that("the error reports how many values are missing", {
  d <- make_data()
  d$Y[c(1, 5, 9)] <- NA
  expect_error(fit_it(d), "Y contains 3 missing values")
})

test_that("missing omics data remains allowed", {
  # The whole point of the incomplete-omics extension: Z may be missing, and
  # every other input may not. A validation that caught Z too would break the
  # package's headline feature, so this is asserted alongside.
  for (regime in c("listwise", "sporadic", "mixed")) {
    d <- make_data(missing = regime)
    expect_silent({
      fit <- fit_it(d)
    })
    expect_s3_class(fit, "early_lucid")
  }
})

test_that("validation applies to parallel and serial too", {
  for (spec in list(list("parallel", 2L, c(2, 2)), list("serial", 2L, c(2, 2)))) {
    model <- spec[[1]]; n_layer <- spec[[2]]; K <- spec[[3]]
    d <- make_data(model, n_layer = n_layer)
    d$Y[3] <- NA
    expect_error(fit_it(d, model = model, K = K), "Y contains 1 missing value",
                 info = model)
  }
})

test_that("validation applies through the lucid() wrapper", {
  d <- make_data()
  d$Y[3] <- NA
  expect_error(
    suppressWarnings(suppressMessages(invisible(capture.output(
      lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
            family = "normal", K = 2, init_omic.data.model = "VVV")
    )))),
    "Y contains 1 missing value"
  )
})

test_that("a binary outcome that is not 0/1 is still rejected on its own terms", {
  # The completeness check runs first, so this must not have started reporting
  # a coding problem as a missingness problem or vice versa.
  d <- make_data()
  d$Y[] <- ifelse(d$Y > 0, 2, 5)
  expect_error(fit_it(d, family = "binary"), "coded as 0 and 1")
})

test_that("initialization arguments are resolved on the fit, not left unevaluated", {
  # estimate_lucid() left init_impute and init_par at their match.arg defaults
  # for the serial path, so a fit taking the defaults reported
  # c("lod", "mix") / c("mclust", "random") -- a record of a choice that was
  # never made, which a bootstrap refit reading these fields would carry forward.
  for (model in c("early", "parallel", "serial")) {
    n_layer <- if (model == "early") 1L else 2L
    K <- if (model == "early") 2 else c(2, 2)
    d <- make_data(model, n_layer = n_layer)
    fit <- fit_it(d, model = model, K = K)

    expect_length(fit$init_impute, 1L)
    expect_length(fit$init_par, 1L)
    expect_equal(fit$init_impute, "lod", info = model)
    expect_equal(fit$init_par, "mclust", info = model)
  }
})

test_that("a requested initialization method is the one recorded", {
  d <- make_data("serial", n_layer = 2L)
  for (ip in c("mix", "lod")) {
    for (pr in c("mclust", "random")) {
      fit <- fit_it(d, model = "serial", K = c(2, 2),
                    init_impute = ip, init_par = pr)
      expect_equal(fit$init_impute, ip)
      expect_equal(fit$init_par, pr)
    }
  }
})

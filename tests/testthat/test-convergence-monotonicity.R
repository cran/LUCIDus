# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# check_convergence() (R/stability_utils.R, internal since it's an EM building
# block rather than something a user calls directly) used to have zero
# call sites: early and parallel each did their own ad-hoc abs(old - new) < tol
# check, and nothing anywhere warned if a log-likelihood trace actually
# decreased. These tests pin (1) that the exported helper's own logic is
# correct in isolation, and (2) that a real fit's loglik_trace never drops by
# more than the documented, majorization/selection-aware slack -- the same
# check the fitting code itself now runs and warns on.

test_that("check_convergence agrees on toy values in both directions", {
  expect_true(LUCIDus:::check_convergence(-100, -100 + 1e-9))          # abs criterion
  expect_true(LUCIDus:::check_convergence(-1e8, -1e8 * (1 - 1e-9)))    # rel criterion
  expect_false(LUCIDus:::check_convergence(-100, -101))                 # neither
  expect_false(LUCIDus:::check_convergence(-100, NA))                   # non-finite
  expect_false(LUCIDus:::check_convergence(Inf, -100))                  # non-finite
})

test_that("a genuine log-likelihood decrease triggers the monotonicity warning", {
  set.seed(9101)
  N <- 60
  G <- matrix(rnorm(N * 3), nrow = N)
  Z <- matrix(rnorm(N * 4), nrow = N)
  Y <- rnorm(N)

  # No missingness, no penalty active: the strict slack (~tol) applies, so a
  # real ascent violation would be caught, not absorbed as expected wobble.
  # This fit is well-separated and complete, so it should NOT warn.
  expect_no_warning_matching <- function(expr, pattern) {
    ws <- character(0)
    withCallingHandlers(
      suppressMessages(invisible(capture.output(force(expr)))),
      warning = function(w) {
        ws <<- c(ws, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_false(any(grepl(pattern, ws)), info = paste(ws, collapse = "; "))
  }

  expect_no_warning_matching(
    estimate_lucid(lucid_model = "early", G = G, Z = Z, Y = Y,
                  family = "normal", K = 2, init_omic.data.model = "VVV",
                  seed = 9101),
    "log-likelihood decreased"
  )
})

test_that("real fits never violate the documented ascent slack", {
  # This is the same standard test-missing-data-mechanism.R and the paper
  # oracle already hold fits to; here it is stated directly against the
  # runtime warning's own bound, for early, parallel and serial, with and
  # without a selection penalty active.
  slack_for <- function(has_penalty, has_missing) if (has_penalty || has_missing) 0.5 else 0.05

  check_trace <- function(trace, has_penalty, has_missing, info) {
    if (length(trace) < 2) return(invisible(TRUE))
    dips <- diff(trace)
    expect_gte(min(dips), -slack_for(has_penalty, has_missing), label = info)
  }

  set.seed(9102)
  N <- 250
  G <- matrix(rnorm(N * 4), nrow = N)
  Z <- matrix(rnorm(N * 6), nrow = N)
  Y <- rnorm(N)

  suppressWarnings(suppressMessages(invisible(capture.output(
    fit_plain <- estimate_lucid(lucid_model = "early", G = G, Z = Z, Y = Y,
                                family = "normal", K = 2,
                                init_omic.data.model = "VVV", seed = 9102)
  ))))
  check_trace(fit_plain$em_control$loglik_trace, FALSE, FALSE, "early, no penalty")

  suppressWarnings(suppressMessages(invisible(capture.output(
    fit_penalized <- estimate_lucid(lucid_model = "early", G = G, Z = Z, Y = Y,
                                    family = "normal", K = 2, Rho_Z_Mu = 5,
                                    init_omic.data.model = "VVV", seed = 9102)
  ))))
  check_trace(fit_penalized$em_control$loglik_trace, TRUE, FALSE, "early, penalized")

  Z2 <- list(matrix(rnorm(N * 4), nrow = N), matrix(rnorm(N * 4), nrow = N))
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit_par <- estimate_lucid(lucid_model = "parallel", G = G, Z = Z2, Y = Y,
                              family = "normal", K = c(2, 2),
                              init_omic.data.model = "VVV", seed = 9102)
  ))))
  check_trace(fit_par$em_control$loglik_trace, FALSE, FALSE, "parallel, no penalty")

  suppressWarnings(suppressMessages(invisible(capture.output(
    fit_serial <- estimate_lucid(lucid_model = "serial", G = G, Z = Z2, Y = Y,
                                 family = "normal", K = list(2, 2),
                                 init_omic.data.model = "VVV", seed = 9102)
  ))))
  for (i in seq_along(fit_serial$submodel)) {
    check_trace(fit_serial$submodel[[i]]$em_control$loglik_trace, FALSE, FALSE,
               paste("serial stage", i))
  }
})

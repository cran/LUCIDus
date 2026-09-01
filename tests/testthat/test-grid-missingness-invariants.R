skip_on_cran()

# Invariants that must hold for each missingness regime, across architectures.
#
# The package branches on a per-subject, per-layer code (R/missing.R:62-66):
#   1 = complete, 2 = sporadic (some features), 3 = listwise (whole layer gone)
# and `impute_flag` is driven by code 2 ALONE. That gives three genuinely
# different regimes, not one parameterised path, and each has invariants the
# others do not:
#
#   sporadic only  -> imputation runs; every row is retained
#   listwise only  -> imputation is skipped entirely; rows are dropped from the
#                     Z-likelihood but still inform G -> X and X -> Y
#   both           -> both, and listwise rows must be re-blanked after the
#                     initial imputation would otherwise have invented values
#
# These are cheap, sharp checks: they fail loudly if the regime dispatch breaks,
# without needing a converged, accurate fit.

N_MISS <- 300L
M_MISS <- 4L

fit_miss <- function(model, K, d, init_impute = "mix", family = "normal") {
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(lucid_model = model, G = d$G, Z = d$Z, Y = d$Y,
                          family = family, K = K, init_impute = init_impute,
                          init_omic.data.model = NULL, seed = 1)
  ))))
  fit
}

test_that("check_na classifies each regime into exactly the expected codes", {
  for (regime in c("none", "listwise", "sporadic", "mixed")) {
    d <- sim_lucid("early", N = N_MISS, K = 2, M = M_MISS,
                   missing = regime, seed = 201)
    na <- check_na(d$Z, lucid_model = "early")
    present <- sort(unique(na$indicator_na))

    expected <- switch(regime,
      none     = 1L,
      listwise = c(1L, 3L),
      sporadic = c(1L, 2L),
      mixed    = c(1L, 2L, 3L))
    expect_equal(present, expected, info = regime)

    # impute_flag tracks SPORADIC rows only -- this is what gates the whole
    # imputation block, so it must not react to listwise missingness.
    expect_equal(na$impute_flag, regime %in% c("sporadic", "mixed"), info = regime)
  }
})

test_that("init_impute has no effect when missingness is purely listwise", {
  # With no sporadic rows, impute_flag is FALSE and the initial-imputation block
  # at R/em.R:170-190 is skipped altogether, so "mix" and "lod" cannot differ.
  # If this ever fails, the regime dispatch has broken and listwise rows are
  # being fed through an imputer.
  d <- sim_lucid("early", N = N_MISS, K = 2, M = M_MISS,
                 missing = "listwise", seed = 202)
  fit_mix <- fit_miss("early", 2, d, init_impute = "mix")
  fit_lod <- fit_miss("early", 2, d, init_impute = "lod")

  expect_equal(fit_mix$likelihood, fit_lod$likelihood)
  expect_equal(unname(fit_mix$res_Mu), unname(fit_lod$res_Mu))
  expect_equal(unname(as.matrix(fit_mix$inclusion.p)),
               unname(as.matrix(fit_lod$inclusion.p)))
})

test_that("init_impute does change the fit when sporadic rows are present", {
  # The mirror of the test above: with sporadic rows the imputer runs, so the
  # two initializers must reach different starting points. A pair of tests that
  # only ever asserted equality would pass just as well on dead code.
  d <- sim_lucid("early", N = N_MISS, K = 2, M = M_MISS,
                 missing = "sporadic", seed = 203)
  fit_mix <- fit_miss("early", 2, d, init_impute = "mix")
  fit_lod <- fit_miss("early", 2, d, init_impute = "lod")

  expect_false(isTRUE(all.equal(fit_mix$Z, fit_lod$Z)))
})

test_that("listwise rows stay NA and observed cells are never altered", {
  for (model in c("early", "parallel", "serial")) {
    n_layer <- if (model == "early") 1L else 2L
    K <- if (model == "early") 2 else rep(2, n_layer)
    for (regime in c("listwise", "sporadic", "mixed")) {
      for (init in c("mix", "lod")) {
        lbl <- paste(model, regime, init)
        d <- sim_lucid(model, N = N_MISS, K = 2, M = M_MISS, n_layer = n_layer,
                       missing = regime, seed = 204)
        fit <- fit_miss(model, K, d, init_impute = init)

        # A serial fit returns the omics as supplied; the imputed data for a
        # stage lives on that stage's sub-model. Check whichever object is the
        # fitted one for this architecture.
        fitted_Z <- if (model == "serial") {
          lapply(fit$submodel, function(s) s$Z)
        } else fit$Z

        expect_observed_unchanged(fitted_Z, d$Z, info = lbl)
        if (regime %in% c("listwise", "mixed")) {
          expect_listwise_preserved(fitted_Z, d$miss_index, info = lbl)
        }
        # sporadic cells, by contrast, must have been filled
        if (regime %in% c("sporadic", "mixed")) {
          zf <- if (is.list(fitted_Z)) fitted_Z[[1]] else fitted_Z
          mi <- if (is.list(d$miss_index)) d$miss_index[[1]] else d$miss_index
          sporadic_rows <- which(rowSums(mi) > 0 & rowSums(mi) < ncol(mi))
          if (length(sporadic_rows)) {
            expect_false(anyNA(zf[sporadic_rows, , drop = FALSE]), info = lbl)
          }
        }
      }
    }
  }
})

test_that("fit$Z is the fitted omics for early and parallel, the input for serial", {
  # This asymmetry is real and easy to trip over: `fit$Z` means "the data the
  # model was fitted to" for early and parallel, but "the data you passed in"
  # for serial, whose imputation happens per stage. Pinned here so a change in
  # either direction is a deliberate one.
  d_e <- sim_lucid("early", N = N_MISS, K = 2, M = M_MISS, missing = "sporadic", seed = 208)
  f_e <- fit_miss("early", 2, d_e)
  expect_true(anyNA(d_e$Z))
  expect_false(anyNA(f_e$Z))

  d_p <- sim_lucid("parallel", N = N_MISS, K = 2, M = M_MISS, n_layer = 2,
                   missing = "sporadic", seed = 208)
  f_p <- fit_miss("parallel", c(2, 2), d_p)
  expect_false(anyNA(f_p$Z[[1]]))

  d_s <- sim_lucid("serial", N = N_MISS, K = 2, M = M_MISS, n_layer = 2,
                   missing = "sporadic", seed = 208)
  f_s <- fit_miss("serial", c(2, 2), d_s)
  expect_true(anyNA(f_s$Z[[1]]))                    # unchanged input
  expect_false(anyNA(f_s$submodel[[1]]$Z))          # imputed, on the sub-model
})

test_that("both missingness mechanisms fit sanely for every architecture", {
  # MCAR vs MAR changes WHICH rows go missing, not the code path, but MAR ties
  # missingness to G and Y and so is the case where a mishandled partition would
  # bias the estimates rather than merely lose precision.
  for (model in c("early", "parallel", "serial")) {
    n_layer <- if (model == "early") 1L else 2L
    K <- if (model == "early") 2 else rep(2, n_layer)
    for (mech in c("MCAR", "MAR")) {
      lbl <- paste(model, mech)
      d <- sim_lucid(model, N = N_MISS, K = 2, M = M_MISS, n_layer = n_layer,
                     missing = "mixed", mechanism = mech, seed = 205)
      fit <- fit_miss(model, K, d)
      expect_sane_fit(fit, info = lbl)
      expect_listwise_preserved(fit$Z, d$miss_index, info = lbl)
    }
  }
})

test_that("a layer missing entirely for some rows does not remove them from G or Y", {
  # BA Eq 11: the Z term drops out of a listwise row's likelihood, but the row
  # still contributes to the exposure and outcome models. The posterior for
  # those rows must therefore still be a proper distribution, driven by G and Y.
  d <- sim_lucid("early", N = N_MISS, K = 2, M = M_MISS,
                 missing = "listwise", miss_ratio = 0.3, seed = 206)
  fit <- fit_miss("early", 2, d)

  lw <- which(rowSums(!d$miss_index) == 0)
  expect_gt(length(lw), 0)
  expect_equal(nrow(fit$inclusion.p), N_MISS)          # no rows dropped
  expect_posterior_valid(fit$inclusion.p[lw, , drop = FALSE])
  # and they are not all pinned to the same cluster, which would indicate the
  # G and Y information was discarded along with the omics
  expect_gt(length(unique(max.col(fit$inclusion.p[lw, , drop = FALSE]))), 1)
})

test_that("per-layer missingness is independent across parallel layers", {
  # One layer being absent for a subject must not blank that subject in another
  # layer: the codes are per (subject, layer), not per subject.
  d <- sim_lucid("parallel", N = N_MISS, K = 2, M = M_MISS, n_layer = 2,
                 missing = "listwise", seed = 207)
  # blank a disjoint set of rows in layer 2
  d$Z[[2]][1:20, ] <- NA
  na <- check_na(d$Z, lucid_model = "parallel")

  expect_length(na$indicator_na, 2L)
  expect_true(all(na$indicator_na[[2]][1:20] == 3L))
  # layer 1's classification for those same rows is decided on its own data
  expect_false(all(na$indicator_na[[1]][1:20] == 3L))
})

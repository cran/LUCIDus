# Heavy: fits multiple LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# estimate_lucid() and tune_lucid() both document `lucid_model` as defaulting
# to "early" (the first element of c("early", "parallel", "serial")). Both
# used to forward that unresolved, length-3 default vector to an internal
# helper (est_lucid()/tune_lucid_auxi()) whose *own* lucid_model default is
# only c("early", "parallel") -- match.arg() there saw a definitely-supplied
# (non-missing), length-3 argument that didn't match its own default and
# rejected it, so calling estimate_lucid()/tune_lucid() without naming
# lucid_model failed outright despite the documented default. Caught by
# actually reading a rendered vignette's output rather than only checking
# that knitting didn't throw.

test_that("estimate_lucid() defaults lucid_model to \"early\" when omitted", {
  set.seed(9001)
  N <- 60
  G <- matrix(rnorm(N * 2), nrow = N)
  Z <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  fit_default <- NULL
  fit_explicit <- NULL
  invisible(capture.output(suppressWarnings({
    fit_default <- estimate_lucid(G = G, Z = Z, Y = Y, family = "normal",
                                  K = 2, max_itr = 3, max_tot.itr = 6, seed = 1)
    fit_explicit <- estimate_lucid(G = G, Z = Z, Y = Y, family = "normal",
                                   lucid_model = "early",
                                   K = 2, max_itr = 3, max_tot.itr = 6, seed = 1)
  })))

  expect_s3_class(fit_default, "early_lucid")
  expect_equal(fit_default$likelihood, fit_explicit$likelihood)
})

test_that("tune_lucid() defaults lucid_model to \"early\" when omitted", {
  set.seed(9002)
  N <- 60
  G <- matrix(rnorm(N * 2), nrow = N)
  Z <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  out <- NULL
  invisible(capture.output(suppressWarnings({
    out <- tune_lucid(G = G, Z = Z, Y = Y, family = "normal", K = 2:3,
                      max_itr = 3, max_tot.itr = 6, seed = 1)
  })))

  expect_true(is.list(out))
  expect_true("tune_list" %in% names(out))
})

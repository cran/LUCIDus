# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# A serial fit runs no EM loop of its own, so its top-level em_control used to
# report only the stopping controls. A user then had no way to see, without
# reaching into every submodel, that a stage had exhausted max_itr. These tests
# pin both the settings (which bootstrap refits reuse) and the aggregated
# diagnostics, and the relationship between them and the submodels.

make_serial_fit <- function(max_itr = 100, seed = 1) {
  set.seed(seed)
  N <- 200
  G <- matrix(rnorm(N * 2), nrow = N)
  colnames(G) <- paste0("g", 1:2)
  Z <- list(matrix(rnorm(N * 4), nrow = N), matrix(rnorm(N * 4), nrow = N))
  Y <- matrix(rnorm(N), ncol = 1)
  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(G = G, Z = Z, Y = Y, lucid_model = "serial",
                          family = "normal", K = list(2, 2),
                          init_omic.data.model = NULL,
                          max_itr = max_itr, seed = seed)
  )))
  fit
}

test_that("serial em_control keeps the stopping controls bootstrap refits reuse", {
  fit <- make_serial_fit()
  expect_true(all(c("tol", "max_itr", "max_tot.itr") %in% names(fit$em_control)))
  expect_equal(fit$em_control$max_itr, 100)
})

test_that("serial em_control reports convergence diagnostics, not just settings", {
  fit <- make_serial_fit()
  ec <- fit$em_control
  expect_true(all(c("converged", "n_iter",
                    "submodel_converged", "submodel_n_iter") %in% names(ec)))
  expect_type(ec$converged, "logical")
  expect_length(ec$converged, 1L)
  expect_length(ec$submodel_converged, length(fit$submodel))
  expect_length(ec$submodel_n_iter, length(fit$submodel))
})

test_that("serial converged is the conjunction over submodels", {
  fit <- make_serial_fit()
  per_stage <- vapply(fit$submodel,
                      function(s) isTRUE(s$em_control$converged), logical(1))
  expect_identical(fit$em_control$converged, all(per_stage))
  expect_equal(unname(fit$em_control$submodel_converged), unname(per_stage))
})

test_that("serial n_iter totals the submodel iteration counts", {
  fit <- make_serial_fit()
  per_stage <- vapply(fit$submodel, function(s) {
    it <- s$em_control$n_iter
    if (is.null(it)) NA_integer_ else as.integer(it)[1]
  }, integer(1))
  expect_equal(fit$em_control$n_iter, sum(per_stage, na.rm = TRUE))
})

test_that("a stage stopped at max_itr makes serial converged FALSE", {
  # max_itr = 1 cannot meet the tolerance, so every stage stops unconverged and
  # the aggregate must say so rather than reporting a converged fit.
  fit <- make_serial_fit(max_itr = 1)
  expect_false(fit$em_control$converged)
  expect_true(any(!fit$em_control$submodel_converged))
})

test_that("serial em_control carries each stage's own loglik_trace", {
  # A serial fit runs no EM loop of its own to trace, so this is a list of
  # each stage's own trace (already computed inside that stage) rather than a
  # single vector the way early/parallel's em_control$loglik_trace is.
  fit <- make_serial_fit()
  ec <- fit$em_control
  expect_true("submodel_loglik_trace" %in% names(ec))
  expect_length(ec$submodel_loglik_trace, length(fit$submodel))
  for (i in seq_along(fit$submodel)) {
    expect_equal(ec$submodel_loglik_trace[[i]],
                fit$submodel[[i]]$em_control$loglik_trace)
    expect_true(is.numeric(ec$submodel_loglik_trace[[i]]))
  }
})

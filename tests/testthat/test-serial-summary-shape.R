# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# The serial summary used to append a second copy of BIC and loglik under the
# same names, leaving the returned list with duplicate entries. `$BIC` resolved
# to the first, so the copies were unreachable, but str() output looked
# malformed and any name-based iteration visited them twice. These tests pin the
# shape, including the legacy alias that downstream code still relies on.

make_serial_fit <- function(seed = 1) {
  set.seed(seed)
  N <- 200
  G <- matrix(rnorm(N * 2), nrow = N)
  colnames(G) <- paste0("g", 1:2)
  Z <- list(matrix(rnorm(N * 4), nrow = N), matrix(rnorm(N * 4), nrow = N))
  Y <- matrix(rnorm(N), ncol = 1)
  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(G = G, Z = Z, Y = Y, lucid_model = "serial",
                          family = "normal", K = list(2, 2),
                          init_omic.data.model = NULL, seed = seed)
  )))
  fit
}

test_that("serial summary has no duplicated component names", {
  s <- summary(make_serial_fit(), auto_print = FALSE)
  expect_false(any(duplicated(names(s))))
})

test_that("serial summary still exposes BIC and loglik at the top level", {
  s <- summary(make_serial_fit(), auto_print = FALSE)
  expect_true(is.numeric(s$BIC) && length(s$BIC) == 1L)
  expect_true(is.numeric(s$loglik) && length(s$loglik) == 1L)
  expect_equal(s$BIC, s$model_fit$BIC)
  expect_equal(s$loglik, s$model_fit$loglik)
})

test_that("summary.list remains a working alias of stage_summary", {
  s <- summary(make_serial_fit(), auto_print = FALSE)
  expect_identical(s$summary.list, s$stage_summary)
  expect_length(s$stage_summary, s$model_info$n_stages)
})

test_that("serial summary still prints without error", {
  fit <- make_serial_fit()
  expect_output(print(summary(fit, auto_print = FALSE)))
})

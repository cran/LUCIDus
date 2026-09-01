# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# The serial results list used to leave `likelihood` and `select` commented
# out (dead references to variables never assigned in that scope) and a `#z =`
# line that was never resurrected. Every submodel already carried its own
# likelihood and selection, so the fix is to surface them at the top level the
# same way early/parallel already do, plus force Rho_G to 0 past stage 1 --
# every stage after the first receives the previous stage's posterior cluster
# probabilities as its "G", not real exposures, so there is nothing there for
# an exposure-selection penalty to select among.

make_penalized_serial_fit <- function(seed = 11, n_layer2_col = 4L) {
  set.seed(seed)
  N <- 300
  G <- matrix(rnorm(N * 6), nrow = N)
  colnames(G) <- paste0("g", 1:6)
  Z1 <- matrix(rnorm(N * 5), nrow = N)
  Z2 <- matrix(rnorm(N * n_layer2_col), nrow = N)
  Y <- rnorm(N)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(G = G, Z = list(Z1, Z2), Y = Y, lucid_model = "serial",
                          family = "normal", K = list(3, 2),
                          init_omic.data.model = "VVV",
                          Rho_G = 0.1, seed = seed)
  ))))
  fit
}

test_that("top-level select is exactly stage 1's own select", {
  fit <- make_penalized_serial_fit()
  expect_identical(fit$select, fit$submodel[[1]]$select)
  expect_true(all(c("selectG", "selectZ") %in% names(fit$select)))
})

test_that("top-level likelihood is the sum of the submodel likelihoods", {
  fit <- make_penalized_serial_fit()
  stage_ll <- vapply(fit$submodel, function(s) s$likelihood, numeric(1))
  expect_true(is.numeric(fit$likelihood))
  expect_length(fit$likelihood, 1L)
  expect_equal(fit$likelihood, sum(stage_ll))
  expect_equal(fit$likelihood, LUCIDus:::cal_loglik_serial(fit))
})

test_that("z is absent from a serial fit", {
  fit <- make_penalized_serial_fit()
  expect_false("z" %in% names(fit))
})

test_that("reorder_lucid and reorder_z no longer exist", {
  expect_false(exists("reorder_lucid", where = asNamespace("LUCIDus"), inherits = FALSE))
  expect_false(exists("reorder_z", where = asNamespace("LUCIDus"), inherits = FALSE))
})

test_that("Rho_G is forced to 0 for every stage after the first, regardless of column count", {
  # Stage 2's incoming "G" is stage 1's non-reference posterior probabilities:
  # K = 3 forwards 2 columns (>= 2, so a column-count guard alone would have
  # left Rho_G nonzero here), while K = 2 forwards exactly 1 (< 2, so even the
  # old buggy guard would have caught this case). Both must show Rho_G == 0.
  fit_k3 <- make_penalized_serial_fit(seed = 21)          # stage 1 K = 3 -> 2 forwarded columns
  expect_equal(fit_k3$submodel[[1]]$Rho$Rho_G, 0.1)
  expect_equal(fit_k3$submodel[[2]]$Rho$Rho_G, 0)

  set.seed(22)
  N <- 300
  G <- matrix(rnorm(N * 6), nrow = N)
  Z1 <- matrix(rnorm(N * 5), nrow = N)
  Z2 <- matrix(rnorm(N * 4), nrow = N)
  Y <- rnorm(N)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit_k2 <- estimate_lucid(G = G, Z = list(Z1, Z2), Y = Y, lucid_model = "serial",
                             family = "normal", K = list(2, 2),
                             init_omic.data.model = "VVV",
                             Rho_G = 0.1, seed = 22)
  ))))
  expect_equal(fit_k2$submodel[[1]]$Rho$Rho_G, 0.1)
  expect_equal(fit_k2$submodel[[2]]$Rho$Rho_G, 0)
})

test_that("a later stage's selectZ can still be a real selection result", {
  # Only selectG is forced meaningless past stage 1; selectZ (omics selection)
  # is unaffected because it is never derived from the pseudo-G.
  fit <- make_penalized_serial_fit()
  expect_type(fit$submodel[[2]]$select$selectZ, "logical")
})

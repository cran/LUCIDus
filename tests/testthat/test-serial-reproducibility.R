# F1.2 / D11: upstream stages of a serial model have no outcome to order by.
# They are now ordered deterministically by their omics cluster means, so the
# reported transition coefficients are reproducible across seeds -- which is
# also what makes bootstrap replicates alignable with the original fit.

test_that("upstream serial stages get a deterministic cluster ordering", {
  skip_on_cran()
  d <- sim_lucid("serial", N = 350, K = c(2, 2), n_layer = 2, M = 3,
                 family = "normal", seed = 81)
  fits <- lapply(c(1, 7, 19, 101), function(s) {
    suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                           K = list(2, 2), family = "normal", seed = s))
  })
  # stage 1 is unsupervised; its clusters must be ordered the same way every run
  ord <- vapply(fits, function(f) {
    mu <- f$submodel[[1]]$res_Mu
    as.integer(mu[1, 1] <= mu[2, 1])
  }, integer(1))
  expect_equal(length(unique(ord)), 1L)

  # and the downstream transition coefficient keeps a stable sign
  dl <- vapply(fits, function(f) f$submodel[[2]]$res_Beta[2, 2], numeric(1))
  expect_equal(length(unique(sign(dl))), 1L)
})

test_that("multiple starts stabilize serial estimates across seeds", {
  skip_on_cran()
  d <- sim_lucid("serial", N = 350, K = c(2, 2), n_layer = 2, M = 3,
                 family = "normal", seed = 81)
  spread <- function(ns) {
    v <- vapply(c(1, 19, 42, 101), function(s) {
      suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                             K = list(2, 2), family = "normal",
                             seed = s, n_starts = ns))$submodel[[2]]$res_Beta[2, 2]
    }, numeric(1))
    diff(range(v))
  }
  s1 <- spread(1); s4 <- spread(4)
  expect_lte(s4, s1 + 1e-8)
  expect_lt(s4, 0.1)
})

test_that("serial summary exposes BIC at the top level like the other models", {
  skip_on_cran()
  d <- sim_lucid("serial", N = 300, K = c(2, 2), n_layer = 2, M = 3,
                 family = "normal", seed = 82)
  f <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                              K = list(2, 2), family = "normal"))
  s <- summary(f, auto_print = FALSE)
  expect_true(is.finite(s$BIC))
  expect_equal(s$BIC, s$model_fit$BIC)
})

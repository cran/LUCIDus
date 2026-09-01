# Variable selection under penalties: does the model keep the exposures and
# omics features that carry signal, and drop the ones that do not?
#
# A penalized fit is deliberately biased, so point recovery is the wrong
# instrument here. These tests assert DIRECTION instead: the right variables
# survive, the null ones are dropped, and estimates shrink toward zero as the
# penalty rises. The truth is known by construction -- sim_default_beta() puts
# the whole exposure effect on the first exposure and leaves the rest at zero,
# and sim_default_mu() separates clusters on every omics feature.

skip_on_cran()

N_SEL <- 800L

fit_pen <- function(d, Rho_G = 0, Rho_Z_Mu = 0, Rho_Z_Cov = 0, K = 2,
                    model = "early", ...) {
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model, family = "normal",
                 K = K, Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu, Rho_Z_Cov = Rho_Z_Cov,
                 init_omic.data.model = "VVV", ...)
  ))))
  fit
}

test_that("Rho_G keeps the signal exposure and drops the null ones", {
  # 1 informative exposure, 5 null ones.
  d <- sim_lucid("early", N = N_SEL, K = 2, M = 4, P = 6, effect = 2, seed = 601)
  fit <- fit_pen(d, Rho_G = 0.05)

  sel <- fit$selection
  expect_false(is.null(sel))
  expect_true(sel$selectG[1], label = "informative exposure retained")
  # the null exposures should mostly go; requiring all six to drop would be
  # brittle, but keeping most of them would mean selection is not working
  expect_lte(sum(sel$selectG), 4L)
})

test_that("raising Rho_G drops more exposures, monotonically", {
  d <- sim_lucid("early", N = N_SEL, K = 2, M = 4, P = 6, effect = 2, seed = 602)
  kept <- vapply(c(0.01, 0.05, 0.2), function(rho) {
    fit <- fit_pen(d, Rho_G = rho)
    if (is.null(fit$selection)) ncol(d$G) else sum(fit$selection$selectG)
  }, numeric(1))

  expect_false(is.unsorted(rev(kept)), label = "kept exposures are non-increasing in Rho_G")
  expect_lt(kept[3], kept[1])
})

test_that("penalized coefficients are shrunk toward zero", {
  # tune_lucid() returns the PENALIZED fit, so its exposure effect must be
  # smaller in magnitude than the unpenalized one. (lucid() refits without a
  # penalty, which is why this reaches for tune_lucid directly.)
  d <- sim_lucid("early", N = N_SEL, K = 2, M = 4, P = 3, effect = 2, seed = 603)
  suppressWarnings(suppressMessages(invisible(capture.output(
    t0 <- tune_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                     family = "normal", K = 2, Rho_G = 0,
                     init_omic.data.model = "VVV")
  ))))
  suppressWarnings(suppressMessages(invisible(capture.output(
    t1 <- tune_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                     family = "normal", K = 2, Rho_G = 0.1,
                     init_omic.data.model = "VVV")
  ))))
  rebase <- function(f) {
    b <- f$best_model$res_Beta
    max(abs(sweep(b, 2, b[1, ], "-")[2, -1]))
  }
  expect_lt(rebase(t1), rebase(t0))
})

test_that("Rho_Z_Mu selects omics features and the informative ones survive", {
  # Half the features separate the clusters, half are pure noise.
  set.seed(604)
  d <- sim_lucid("early", N = N_SEL, K = 2, M = 4, seed = 604)
  noise <- matrix(rnorm(N_SEL * 4), nrow = N_SEL)
  colnames(noise) <- paste0("noise", 1:4)
  d$Z <- cbind(d$Z, noise)

  fit <- fit_pen(d, Rho_Z_Mu = 20)
  sel <- fit$selection
  expect_false(is.null(sel))
  # the four generating features carry a mean difference of 2*sep; the four
  # noise columns carry none, so selection must favour the former
  expect_gte(sum(sel$selectZ[1:4]), sum(sel$selectZ[5:8]))
})

test_that("Rho_Z_Cov alone engages the graphical-lasso path", {
  # Rho_Z_Mu = 0 with Rho_Z_Cov > 0 is a distinct branch: penalized covariance
  # estimation runs while mean soft-thresholding does not.
  d <- sim_lucid("early", N = N_SEL, K = 2, M = 4, seed = 605)
  fit <- fit_pen(d, Rho_Z_Mu = 0, Rho_Z_Cov = 0.1)

  expect_sane_fit(fit, info = "Rho_Z_Cov only")
  expect_equal(fit$Rho$Rho_Z_Cov, 0.1)
  expect_equal(fit$Rho$Rho_Z_Mu, 0)
  # clusters must still be recovered despite the covariance penalty
  expect_gt(cluster_accuracy(fit$inclusion.p, d$truth$X[[1]]), 0.85)
})

test_that("an unpenalized fit selects everything", {
  d <- sim_lucid("early", N = N_SEL, K = 2, M = 4, P = 3, seed = 606)
  fit <- fit_pen(d, Rho_G = 0, Rho_Z_Mu = 0, Rho_Z_Cov = 0)
  expect_true(all(fit$select$selectG))
  expect_true(all(fit$select$selectZ))
})

test_that("lucid() refits the selected model without a penalty", {
  # The refit is what makes lucid()'s estimates unshrunk, and it is the reason
  # a `selection` component exists separately from `select`.
  d <- sim_lucid("early", N = N_SEL, K = 2, M = 4, P = 6, effect = 2, seed = 607)
  fit <- fit_pen(d, Rho_G = 0.05)

  expect_false(is.null(fit$selection))
  # `select` describes the refit dimensions, `selection` the original inputs
  expect_equal(length(fit$selection$selectG), ncol(d$G))
  expect_equal(ncol(fit$res_Beta) - 1L, sum(fit$selection$selectG))
  # `Rho` keeps the tuned penalty as metadata -- it records what produced the
  # selection, not what the refit ran with, and `selection$Rho` agrees.
  expect_equal(fit$Rho$Rho_G, 0.05)
  expect_equal(fit$selection$Rho$Rho_G, 0.05)

  # What makes the refit unpenalized is that its estimates are not shrunk: on
  # the retained exposure the effect must be at least as large as the penalized
  # fit produced on the same data.
  suppressWarnings(suppressMessages(invisible(capture.output(
    tuned <- tune_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                        family = "normal", K = 2, Rho_G = 0.05,
                        init_omic.data.model = "VVV")
  ))))
  pen_b <- tuned$best_model$res_Beta
  pen_effect <- max(abs(sweep(pen_b, 2, pen_b[1, ], "-")[2, -1]))
  ref_b <- fit$res_Beta
  ref_effect <- max(abs(sweep(ref_b, 2, ref_b[1, ], "-")[2, -1]))
  expect_gte(ref_effect, pen_effect)
})

test_that("parallel selection reports per-layer and union exposure sets", {
  d <- sim_lucid("parallel", N = N_SEL, K = 2, M = 4, P = 4, n_layer = 2,
                 effect = 2, seed = 608)
  fit <- fit_pen(d, Rho_G = 0.05, K = list(2, 2), model = "parallel")

  expect_false(is.null(fit$select$selectG_layer))
  expect_length(fit$select$selectG_layer, 2L)
  # the overall set is the union across layers: anything selected in a layer
  # must appear in it
  union_layers <- Reduce(`|`, fit$select$selectG_layer)
  expect_equal(as.logical(fit$select$selectG), as.logical(union_layers))
})

test_that("penalized fits still recover clusters, if less precisely", {
  # Graceful degradation: shrinkage costs accuracy in the coefficients but must
  # not destroy the clustering the model exists to find.
  d <- sim_lucid("early", N = N_SEL, K = 2, M = 4, P = 4, effect = 2, seed = 609)
  for (rho in c(0, 0.05, 0.15)) {
    fit <- fit_pen(d, Rho_G = rho)
    expect_gt(cluster_accuracy(fit$inclusion.p, d$truth$X[[1]]), 0.85,
              label = paste("cluster accuracy at Rho_G =", rho))
  }
})

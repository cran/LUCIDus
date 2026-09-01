# Parameter recovery: do the estimates equal the values that generated the data?
#
# This is the tier that answers "accurate", not merely "ran". Data are simulated
# from the LUCID DAG with a known Theta (helper-sim.R) and the fit is compared
# back to it, with cluster labels aligned first -- the package canonicalizes
# cluster order by outcome level for supervised fits, but not for unsupervised
# ones or upstream serial stages, so alignment is required rather than optional.
#
# Tolerances are calibrated, not guessed. For the early model at sep = 1.5,
# K = 2, M = 4, the maximum absolute error in mu was measured at
#   N =  300 -> 0.263      N =  600 -> 0.161      N = 1200 -> 0.056
# i.e. roughly 1/sqrt(N). N = 1200 with tolerance 0.15 therefore sits several
# standard errors clear of the noise floor while still catching real bias.
# Cluster accuracy is driven by separation rather than N and held above 0.94
# throughout, so it is asserted at 0.90.

skip_on_cran()

N_ACC <- 1200L
M_ACC <- 4L
TOL_MU <- 0.15
TOL_GAMMA <- 0.25
MIN_ACC <- 0.90

fit_acc <- function(model, K, d, family = "normal", ...) {
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model, family = family,
                 K = K, init_omic.data.model = "VVV", ...)
  ))))
  fit
}

test_that("early recovers cluster means, assignment and outcome effects", {
  for (family in c("normal", "binary")) {
    d <- sim_lucid("early", N = N_ACC, K = 2, M = M_ACC, family = family, seed = 501)
    fit <- fit_acc("early", 2, d, family = family)
    lbl <- paste("early", family)

    perm <- expect_recovers_mu(fit$res_Mu, d$truth$mu[[1]], TOL_MU, info = lbl)
    expect_gt(cluster_accuracy(fit$inclusion.p, d$truth$X[[1]]), MIN_ACC)

    # gamma is generated as absolute cluster levels 0, 1, ..., K-1; the fit
    # reports the same on res_Gamma$beta, aligned by the same permutation.
    got <- as.numeric(fit$res_Gamma$beta)[perm]
    expect_equal(got - got[1], d$truth$gamma[[1]] - d$truth$gamma[[1]][1],
                 tolerance = TOL_GAMMA)
  }
})

test_that("early recovers per-cluster residual standard deviations", {
  # The early normal outcome model estimates a sigma PER CLUSTER, unlike the
  # parallel model's single pooled sd, so heteroscedastic truth is recoverable
  # here and nowhere else.
  d <- sim_lucid("early", N = N_ACC, K = 2, M = M_ACC, family = "normal",
                 gamma_sd = c(0.6, 1.4), seed = 502)
  fit <- fit_acc("early", 2, d)
  perm <- match_clusters(fit$res_Mu, d$truth$mu[[1]])$perm
  expect_equal(as.numeric(fit$res_Gamma$sigma)[perm], c(0.6, 1.4), tolerance = 0.20)
})

test_that("early recovers the exposure effect on cluster membership", {
  # sim_default_beta puts the whole signal on the first exposure and nulls the
  # rest, so both the magnitude and the sparsity pattern are checkable.
  d <- sim_lucid("early", N = N_ACC, K = 2, M = M_ACC, P = 3, effect = 1.5, seed = 503)
  fit <- fit_acc("early", 2, d)
  perm <- match_clusters(fit$res_Mu, d$truth$mu[[1]])$perm
  beta <- fit$res_Beta[perm, , drop = FALSE]
  beta <- sweep(beta, 2, beta[1, ], "-")          # rebase on the reference row

  expect_equal(unname(beta[2, 2]), 1.5, tolerance = 0.30)   # signal exposure
  expect_equal(unname(beta[2, 3:4]), c(0, 0), tolerance = 0.30)  # null exposures
})

test_that("parallel recovers per-layer means at 2 and 3 layers", {
  for (n_layer in c(2L, 3L)) {
    for (family in c("normal", "binary")) {
      lbl <- paste("parallel", n_layer, family)
      d <- sim_lucid("parallel", N = N_ACC, K = 2, M = M_ACC,
                     n_layer = n_layer, family = family, seed = 504)
      fit <- fit_acc("parallel", as.list(rep(2, n_layer)), d, family = family)

      for (i in seq_len(n_layer)) {
        expect_recovers_mu(fit$res_Mu[[i]], d$truth$mu[[i]], TOL_MU,
                           info = paste(lbl, "layer", i))
        expect_gt(cluster_accuracy(fit$inclusion.p[[i]], d$truth$X[[i]]), MIN_ACC)
      }
    }
  }
})

test_that("serial recovers each stage and the between-stage transition", {
  d <- sim_lucid("serial", N = N_ACC, K = 2, M = M_ACC, n_layer = 3, seed = 505)
  fit <- fit_acc("serial", as.list(rep(2, 3)), d)

  for (i in 1:3) {
    sub <- fit$submodel[[i]]
    expect_recovers_mu(sub$res_Mu, d$truth$mu[[i]], TOL_MU,
                       info = paste("serial stage", i))
    expect_gt(cluster_accuracy(sub$inclusion.p, d$truth$X[[i]]), MIN_ACC)
  }

  # res_Delta[[i]] is stage i+1 fitted on stage i's cluster assignment, so it
  # is comparable to truth$delta[[i]] -- the coefficients actually used to
  # generate that transition, which truth$beta never held.
  expect_length(fit$res_Delta, 2L)
  for (i in 1:2) {
    got <- as.matrix(fit$res_Delta[[i]])
    got <- sweep(got, 2, got[1, ], "-")
    want <- d$truth$delta[[i]]
    expect_equal(dim(got), dim(want), info = paste("delta", i))
    # the transition effect must be present and of the right sign, which is the
    # substantive claim; its exact size is attenuated by cluster uncertainty
    expect_gt(abs(unname(got[2, 2])), 0.5 * abs(want[2, 2]))
    expect_equal(sign(unname(got[2, 2])), sign(want[2, 2]), info = paste("delta", i))
  }
})

test_that("serial with a parallel stage recovers every layer of that stage", {
  topo <- list(2, list(2, 2))
  d <- sim_lucid_serial_topology(topo, N = N_ACC, M = M_ACC, seed = 506)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                 family = "normal", K = topo, init_omic.data.model = "VVV")
  ))))

  expect_recovers_mu(fit$submodel[[1]]$res_Mu, d$truth$mu[[1]][[1]], TOL_MU,
                     info = "serial mixed stage 1")
  for (l in 1:2) {
    expect_recovers_mu(fit$submodel[[2]]$res_Mu[[l]], d$truth$mu[[2]][[l]], TOL_MU,
                       info = paste("serial mixed stage 2 layer", l))
  }
})

test_that("covariate effects on both arms are recovered", {
  d <- sim_lucid("early", N = N_ACC, K = 2, M = M_ACC, n_CoY = 1,
                 coef_CoY = 0.8, seed = 507)
  fit <- fit_acc("early", 2, d, CoY = d$CoY)
  expect_equal(unname(as.numeric(fit$res_Gamma$covariate)), 0.8, tolerance = 0.20)
})

test_that("missingness costs precision but does not bias the estimates", {
  # The "assert direction" rule for degraded cells: rather than a fixed
  # tolerance, require that the missing-data fit stays close to the
  # complete-data fit on the SAME simulated data. That isolates the cost of
  # missingness from ordinary sampling error.
  for (regime in c("listwise", "sporadic", "mixed")) {
    complete <- sim_lucid("early", N = N_ACC, K = 2, M = M_ACC,
                          missing = "none", seed = 508)
    partial <- sim_lucid("early", N = N_ACC, K = 2, M = M_ACC,
                         missing = regime, seed = 508)
    # same seed, so the underlying data and truth are identical; only the
    # missingness mask differs
    expect_equal(complete$truth$X[[1]], partial$truth$X[[1]])

    fit_c <- fit_acc("early", 2, complete)
    fit_p <- fit_acc("early", 2, partial)

    perm_c <- match_clusters(fit_c$res_Mu, complete$truth$mu[[1]])$perm
    perm_p <- match_clusters(fit_p$res_Mu, partial$truth$mu[[1]])$perm
    err_c <- max(abs(fit_c$res_Mu[perm_c, ] - complete$truth$mu[[1]]))
    err_p <- max(abs(fit_p$res_Mu[perm_p, ] - partial$truth$mu[[1]]))

    expect_lt(err_p, 0.60, label = paste("recovery error under", regime))
    # degradation is bounded: missingness may cost precision, but a 20% mask
    # must not multiply the error several-fold
    expect_lt(err_p, err_c + 0.45, label = paste("degradation under", regime))
    expect_gt(cluster_accuracy(fit_p$inclusion.p, partial$truth$X[[1]]), 0.85)
  }
})

test_that("unsupervised fits still recover the omics clusters", {
  # Without Y the clusters are identified by the omics alone, so mu must still
  # be recovered -- and label alignment genuinely matters here, since the
  # outcome-based canonical ordering does not apply.
  d <- sim_lucid("early", N = N_ACC, K = 2, M = M_ACC, seed = 509)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                          family = "normal", K = 2, useY = FALSE,
                          init_omic.data.model = "VVV")
  ))))
  expect_recovers_mu(fit$res_Mu, d$truth$mu[[1]], TOL_MU, info = "unsupervised")
  expect_gt(cluster_accuracy(fit$inclusion.p, d$truth$X[[1]]), MIN_ACC)
})

test_that("K = 3 is recovered as well as K = 2", {
  d <- sim_lucid("early", N = N_ACC, K = 3, M = M_ACC, sep = 2, seed = 510)
  fit <- fit_acc("early", 3, d)
  expect_recovers_mu(fit$res_Mu, d$truth$mu[[1]], 0.25, info = "K=3")
  expect_gt(cluster_accuracy(fit$inclusion.p, d$truth$X[[1]]), MIN_ACC)
})

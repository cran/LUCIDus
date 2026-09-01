# Heavy: fits multiple LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# F1: relabelling invariance, and "every reported block refers to the same
# cluster".  The 3.0.3 gamma display bug -- correct estimate, wrong label -- is
# exactly the failure class these tests exist to catch.

test_that("the likelihood is invariant to a cluster relabelling", {
  d <- sim_lucid("early", N = 400, K = 2, M = 4, family = "normal", seed = 51)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  S <- simplify2array(fit$res_Sigma)
  base <- o_estep_early(d$G, d$Z, d$Y, fit$res_Beta, fit$res_Mu, S,
                        fit$res_Gamma$cluster_effect, fit$res_Gamma$sigma, "normal", TRUE)
  perm <- c(2, 1)
  bperm <- fit$res_Beta[perm, , drop = FALSE]
  bperm <- sweep(bperm, 2, bperm[1, ], "-")
  sw <- o_estep_early(d$G, d$Z, d$Y, bperm, fit$res_Mu[perm, , drop = FALSE], S[, , perm],
                      fit$res_Gamma$cluster_effect[perm], fit$res_Gamma$sigma[perm],
                      "normal", TRUE)
  expect_equal(base$loglik, sw$loglik, tolerance = 1e-8)
  expect_equal(unname(base$r), unname(sw$r[, perm]), tolerance = 1e-10)
})

test_that("clusters are ordered by outcome level and beta is rebased on cluster 1", {
  for (fam in c("normal", "binary")) {
    d <- sim_lucid("early", N = 400, K = 3, M = 4, family = fam, seed = 52)
    fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                  family = fam, K = 3, init_omic.data.model = "VVV"))
    expect_false(is.unsorted(fit$res_Gamma$cluster_effect), info = fam)
    expect_true(all(fit$res_Beta[1, ] == 0), info = fam)
  }
})

test_that("gamma contrasts reconstruct the absolute cluster levels", {
  # This is the 3.0.3 regression: beta held absolute means but printed as if
  # reference-coded, so the reported cluster-2 "effect" was its own mean.
  for (K in c(2, 3)) {
    d <- sim_lucid("early", N = 400, K = K, M = 4, family = "normal", seed = 53)
    fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                  family = "normal", K = K, init_omic.data.model = "VVV"))
    b <- as.numeric(fit$res_Gamma$beta)
    rebuilt <- c(b[1], b[1] + b[seq.int(2, K)])
    expect_equal(rebuilt, as.numeric(fit$res_Gamma$cluster_effect), tolerance = 1e-8)
    # the first coefficient is cluster 1's LEVEL, not zero
    expect_false(isTRUE(all.equal(b[1], 0)))
  }
})

test_that("inclusion.p column k corresponds to mu row k", {
  d <- sim_lucid("early", N = 500, K = 3, M = 4, family = "normal", seed = 54)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 3, init_omic.data.model = "VVV"))
  r <- fit$inclusion.p
  expect_equal(unname((t(r) %*% d$Z) / colSums(r)), unname(fit$res_Mu), tolerance = 1e-2)
})

test_that("summary tables move together with the cluster labels", {
  d <- sim_lucid("early", N = 400, K = 2, M = 4, family = "normal", seed = 55)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  s <- summary(fit, auto_print = FALSE)
  # mu table column k is res_Mu row k
  expect_equal(unname(as.matrix(s$parameters$mu)), unname(fit$res_Mu), tolerance = 1e-12)
  # gamma table is the reference-coded vector
  expect_equal(as.numeric(s$parameters$gamma$beta), as.numeric(fit$res_Gamma$beta))
})

test_that("parallel relabelling permutes Beta, Mu, Sigma and the posterior together", {
  d <- sim_lucid("parallel", N = 300, K = c(2, 2), n_layer = 2, M = 3,
                 family = "normal", seed = 56)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, K = list(2, 2),
                                lucid_model = "parallel", family = "normal"))
  for (i in 1:2) {
    r <- fit$inclusion.p[[i]]
    expect_posterior_valid(r)
    # weighted means of the layer's own Z recover that layer's Mu
    expect_equal(unname((t(r) %*% fit$Z[[i]]) / colSums(r)),
                 unname(fit$res_Mu[[i]]), tolerance = 5e-2, info = paste("layer", i))
    expect_equal(nrow(fit$res_Beta$Beta[[i]]), fit$K[i] - 1L)
  }
})

test_that("multi-start reports the best optimum, not an arbitrary one (D10)", {
  skip_on_cran()
  d <- sim_lucid("early", N = 400, K = 2, M = 4, family = "normal", seed = 57)
  # Multi-start's first start reuses the supplied seed, so it can never do worse
  # than the corresponding single-start fit.  (It is not comparable to an
  # arbitrary OTHER seed, which may happen to find a better optimum.)
  single <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                   family = "normal", K = 2,
                                   init_omic.data.model = "VVV", seed = 123))
  multi <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                  family = "normal", K = 2, init_omic.data.model = "VVV",
                                  seed = 123, n_starts = 5))
  expect_gte(multi$likelihood, single$likelihood - 1e-6)
  expect_equal(multi$likelihood, max(multi$em_control$start_loglik, na.rm = TRUE))
  expect_equal(multi$em_control$n_starts, 5L)
  expect_true(multi$em_control$n_starts_ok >= 1L)

  # and it must actually be selecting: with a spread across starts, the chosen
  # fit is strictly better than the worst one
  sl <- multi$em_control$start_loglik
  skip_if(diff(range(sl, na.rm = TRUE)) < 1e-6, "starts all reached the same optimum")
  expect_gt(multi$likelihood, min(sl, na.rm = TRUE))
})

test_that("convergence status is reported honestly (D10)", {
  d <- sim_lucid("early", N = 300, K = 2, M = 4, family = "normal", seed = 58)
  ok <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                               family = "normal", K = 2, init_omic.data.model = "VVV"))
  expect_true(ok$em_control$converged)
  expect_gt(ok$em_control$n_iter, 0)

  expect_warning(
    bad <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                 family = "normal", K = 2, init_omic.data.model = "VVV",
                 max_itr = 2, max_tot.itr = 4),
    "max_itr"
  )
  expect_false(bad$em_control$converged)
})

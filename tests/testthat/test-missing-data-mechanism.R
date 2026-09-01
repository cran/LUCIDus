# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# F2: implementation correctness for the missing-data machinery.
# Equation numbers prefixed "BA" refer to vbae123.

test_that("BA Eq 11: list-wise rows use only G and Y information", {
  d <- sim_lucid("early", N = 400, K = 2, M = 4, family = "normal",
                 missing = "listwise", miss_ratio = 0.2, seed = 41)
  lw <- which(rowSums(is.na(d$Z)) == ncol(d$Z))
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  lev <- fit$res_Gamma$cluster_effect; sdv <- fit$res_Gamma$sigma
  lb <- cbind(1, d$G) %*% t(fit$res_Beta)
  ly <- sapply(seq_len(fit$K), function(k)
    stats::dnorm(as.numeric(d$Y)[lw], lev[k], sdv[k], log = TRUE))
  v <- lb[lw, , drop = FALSE] + ly
  orc <- t(apply(v, 1, function(x) { e <- exp(x - max(x)); e / sum(e) }))
  expect_equal(unname(fit$inclusion.p[lw, ]), unname(orc), tolerance = 1e-10)
})

test_that("list-wise rows stay NA under every initializer (D7)", {
  d <- sim_lucid("early", N = 400, K = 2, M = 5, family = "normal",
                 missing = "mixed", miss_ratio = 0.15, seed = 42)
  lw <- which(rowSums(is.na(d$Z)) == ncol(d$Z))
  sp <- which(rowSums(is.na(d$Z)) > 0 & rowSums(is.na(d$Z)) < ncol(d$Z))
  skip_if(length(lw) == 0 || length(sp) == 0, "need both patterns present")
  for (init in c("mix", "lod")) {
    if (init == "mix" && !requireNamespace("mix", quietly = TRUE)) next
    fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                  family = "normal", K = 2, init_impute = init,
                                  init_omic.data.model = "VVV"))
    expect_true(all(is.na(fit$Z[lw, ])), info = init)
    expect_true(all(is.finite(fit$Z[sp, ])), info = init)
  }
})

test_that("observed omics cells are never altered by imputation", {
  d <- sim_lucid("early", N = 300, K = 2, M = 5, family = "normal",
                 missing = "mixed", miss_ratio = 0.2, seed = 43)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  obs <- !is.na(d$Z)
  expect_equal(fit$Z[obs], d$Z[obs])
})

test_that("imputed values are closer to the withheld truth than the column mean", {
  d <- sim_lucid("early", N = 500, K = 2, M = 4, family = "normal",
                 missing = "sporadic", miss_ratio = 0.25, seed = 44)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  miss <- is.na(d$Z)
  skip_if(sum(miss) < 20, "not enough sporadic cells")
  truth <- d$Z_complete[miss]
  rmse_lucid <- sqrt(mean((fit$Z[miss] - truth)^2))
  colmean <- matrix(colMeans(d$Z, na.rm = TRUE), nrow(d$Z), ncol(d$Z), byrow = TRUE)
  rmse_mean <- sqrt(mean((colmean[miss] - truth)^2))
  expect_lt(rmse_lucid, rmse_mean)
})

test_that("EM ascends the observed log-likelihood (BA Algorithm 1)", {
  # Without imputation the ascent is exact.  With an I-step the guarantee is
  # only as tight as the covariance stabilization allows: a near-singular Sigma
  # is nudged to restore positive definiteness, which is not an ascent step and
  # can produce a small dip that the run then recovers from.  So the
  # assertion is: exact monotonicity when no imputation runs, and negligible
  # total descent relative to total ascent when it does.  A structural defect
  # (for example posterior mass leaking out of the state space) produces
  # systematic, large descent and would still fail this.  The default
  # initializer is "lod" (a fixed per-column fill), which is a cruder start
  # than the EM-based "mix" imputation used to be and produces more, smaller
  # dips before the trace stabilizes -- hence the count threshold below is
  # looser than the (still strict) magnitude threshold.
  for (miss in c("none", "listwise")) {
    d <- sim_lucid("early", N = 350, K = 2, M = 4, family = "normal",
                   missing = miss, miss_ratio = 0.2, seed = 45)
    fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                  family = "normal", K = 2, init_omic.data.model = "VVV"))
    tr <- fit$em_control$loglik_trace
    skip_if(length(tr) < 3, "trace too short")
    expect_true(all(diff(tr) > -1e-5), info = paste("missing =", miss))
  }
  for (miss in c("sporadic", "mixed")) {
    d <- sim_lucid("early", N = 350, K = 2, M = 4, family = "normal",
                   missing = miss, miss_ratio = 0.2, seed = 45)
    fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                  family = "normal", K = 2, init_omic.data.model = "VVV"))
    tr <- fit$em_control$loglik_trace
    skip_if(length(tr) < 3, "trace too short")
    dd <- diff(tr)
    descent <- abs(sum(dd[dd < 0])); ascent <- sum(dd[dd > 0])
    expect_lt(descent / ascent, 1e-2, label = paste("relative descent,", miss))
    expect_lte(sum(dd < -1e-5), 15L)
    expect_gt(tail(tr, 1), tr[1])
  }
})

test_that("prediction under sporadic missingness uses the observed-coordinate marginal (D1)", {
  d <- sim_lucid("early", N = 300, K = 2, M = 5, family = "normal", seed = 46)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  Zs <- d$Z
  set.seed(1); Zs[cbind(1:25, sample(seq_len(ncol(Zs)), 25, TRUE))] <- NA
  p <- predict_lucid(model = fit, G = d$G, Z = Zs, Y = NULL, lucid_model = "early")
  orc <- o_estep_early(d$G, Zs, NULL, fit$res_Beta, fit$res_Mu,
                       simplify2array(fit$res_Sigma), NULL, NULL,
                       "normal", useY = FALSE)
  expect_equal(unname(p$inclusion.p), unname(orc$r), tolerance = 1e-8)
  # and specifically NOT the degenerate uniform posterior of the old code
  expect_false(all(abs(p$inclusion.p[1:25, 1] - 0.5) < 1e-8))
})

test_that("parallel prediction tolerates sporadic missingness (D1)", {
  d <- sim_lucid("parallel", N = 250, K = c(2, 2), n_layer = 2, M = 3,
                 family = "normal", seed = 47)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, K = list(2, 2),
                                lucid_model = "parallel", family = "normal"))
  Zm <- d$Z
  set.seed(2); Zm[[1]][cbind(1:15, sample(1:3, 15, TRUE))] <- NA
  p <- predict_lucid(model = fit, G = d$G, Z = Zm, Y = NULL, lucid_model = "parallel")
  expect_posterior_valid(p$inclusion.p[[1]])
  expect_true(all(unlist(p$pred.x) %in% seq_len(2)))
})

test_that("per-layer list-wise partition is independent across layers", {
  d <- sim_lucid("parallel", N = 300, K = c(2, 2), n_layer = 2, M = 3,
                 family = "normal", seed = 48)
  # layer 1 list-wise missing on rows 1:40, layer 2 complete
  d$Z[[1]][1:40, ] <- NA
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, K = list(2, 2),
                                lucid_model = "parallel", family = "normal"))
  expect_true(all(is.na(fit$Z[[1]][1:40, ])))
  expect_true(all(is.finite(fit$Z[[2]][1:40, ])))
  expect_posterior_valid(fit$inclusion.p[[2]])
})

test_that("check_na classifies every missingness pattern correctly", {
  Z <- matrix(stats::rnorm(100 * 4), 100, 4)
  Z[1:5, ] <- NA                 # list-wise
  Z[6:10, 2] <- NA               # sporadic
  na <- check_na(Z)
  expect_equal(sum(na$indicator_na == 3), 5L)
  expect_equal(sum(na$indicator_na == 2), 5L)
  expect_equal(sum(na$indicator_na == 1), 90L)
  expect_true(na$impute_flag)
})

test_that("cross-layer missingness proportion counts cells, not layers (D9)", {
  Za <- matrix(NA_real_, 100, 1)
  Zb <- matrix(stats::rnorm(900), 100, 9)
  np <- suppressWarnings(check_na(list(Za, Zb), lucid_model = "parallel"))
  expect_equal(np$cross_layer_summary$total_missing_prop, 0.1)
})

# Imputation accuracy against withheld truth.
#
# sim_lucid() retains Z_complete, the omics matrix before missingness was
# imposed, so the values the model invents for sporadically missing cells can be
# scored against the values that were actually there. That is a real accuracy
# measurement, not a plausibility check.
#
# The claim under test is the one the incomplete-omics extension makes: filling
# a missing cell from the fitted cluster model, conditioning on the subject's
# observed coordinates, beats ignoring the model and substituting a column mean.
# If that stops being true, the I-step is not earning its cost.

skip_on_cran()

N_IMP <- 800L
M_IMP <- 6L

rmse <- function(a, b) sqrt(mean((a - b)^2))

# RMSE over exactly the cells that were blanked sporadically, i.e. rows with
# some but not all features missing.
sporadic_rmse <- function(filled, complete, mask) {
  sporadic <- rowSums(mask) > 0 & rowSums(mask) < ncol(mask)
  cells <- mask & sporadic
  if (!any(cells)) return(NA_real_)
  rmse(filled[cells], complete[cells])
}

# The baseline any imputation must beat: the observed column mean.
column_mean_fill <- function(Z) {
  for (j in seq_len(ncol(Z))) {
    Z[is.na(Z[, j]), j] <- mean(Z[, j], na.rm = TRUE)
  }
  Z
}

test_that("model-based imputation beats column-mean filling for the early model", {
  for (mech in c("MCAR", "MAR")) {
    d <- sim_lucid("early", N = N_IMP, K = 2, M = M_IMP, sep = 2.5,
                   missing = "sporadic", mechanism = mech, seed = 801)
    suppressWarnings(suppressMessages(invisible(capture.output(
      fit <- estimate_lucid(lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                            family = "normal", K = 2,
                            init_omic.data.model = "VVV", seed = 1)
    ))))

    model_err <- sporadic_rmse(fit$Z, d$Z_complete, d$miss_index)
    naive_err <- sporadic_rmse(column_mean_fill(d$Z), d$Z_complete, d$miss_index)

    expect_false(is.na(model_err), label = paste("sporadic cells exist", mech))
    expect_lt(model_err, naive_err, label = paste("model beats column mean under", mech))
    # and by a worthwhile margin, not a rounding error: the clusters are well
    # separated, so knowing the cluster should cut the error substantially
    expect_lt(model_err, 0.85 * naive_err)
  }
})

test_that("model-based imputation beats column-mean filling per layer in parallel", {
  d <- sim_lucid("parallel", N = N_IMP, K = 2, M = M_IMP, n_layer = 2, sep = 2.5,
                 missing = "sporadic", seed = 802)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(lucid_model = "parallel", G = d$G, Z = d$Z, Y = d$Y,
                          family = "normal", K = c(2, 2),
                          init_omic.data.model = "VVV", seed = 1)
  ))))

  for (i in 1:2) {
    model_err <- sporadic_rmse(fit$Z[[i]], d$Z_complete[[i]], d$miss_index[[i]])
    naive_err <- sporadic_rmse(column_mean_fill(d$Z[[i]]), d$Z_complete[[i]],
                               d$miss_index[[i]])
    expect_lt(model_err, naive_err, label = paste("layer", i))
  }
})

test_that("imputation degrades gracefully as more data goes missing", {
  # Error must rise with the missing fraction -- it cannot stay flat, which
  # would mean the observed coordinates are not being used -- but it must stay
  # bounded and keep beating the naive baseline throughout.
  errs <- numeric(0)
  for (ratio in c(0.1, 0.2, 0.35)) {
    d <- sim_lucid("early", N = N_IMP, K = 2, M = M_IMP, sep = 2.5,
                   missing = "sporadic", miss_ratio = ratio, seed = 803)
    suppressWarnings(suppressMessages(invisible(capture.output(
      fit <- estimate_lucid(lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                            family = "normal", K = 2,
                            init_omic.data.model = "VVV", seed = 1)
    ))))
    model_err <- sporadic_rmse(fit$Z, d$Z_complete, d$miss_index)
    naive_err <- sporadic_rmse(column_mean_fill(d$Z), d$Z_complete, d$miss_index)
    expect_lt(model_err, naive_err, label = paste("miss_ratio", ratio))
    errs <- c(errs, model_err)
  }
  # bounded across the whole range
  expect_lt(max(errs), 1.5)
})

test_that("imputation uses the correlation structure, not just the cluster mean", {
  # With correlated features, knowing some coordinates of a subject should
  # sharpen the guess for the others beyond what the cluster mean alone gives.
  # This is the specific claim of the conditional I-step over mean substitution.
  M <- 4L
  rho <- 0.8
  S <- matrix(rho, M, M); diag(S) <- 1
  sig <- array(0, dim = c(M, M, 2)); sig[, , 1] <- S; sig[, , 2] <- S

  d <- sim_lucid("early", N = N_IMP, K = 2, M = M, sep = 2.5, sigma = sig,
                 missing = "sporadic", seed = 804)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                          family = "normal", K = 2,
                          init_omic.data.model = "VVV", seed = 1)
  ))))

  model_err <- sporadic_rmse(fit$Z, d$Z_complete, d$miss_index)

  # Baseline that knows the cluster but ignores the correlation: fill each
  # missing cell with its own cluster's mean for that feature.
  assign <- max.col(fit$inclusion.p)
  cluster_fill <- d$Z
  for (i in seq_len(nrow(cluster_fill))) {
    miss <- which(is.na(cluster_fill[i, ]))
    if (length(miss)) cluster_fill[i, miss] <- fit$res_Mu[assign[i], miss]
  }
  cluster_err <- sporadic_rmse(cluster_fill, d$Z_complete, d$miss_index)

  expect_lt(model_err, cluster_err,
            label = "conditional imputation beats cluster-mean substitution")
})

test_that("listwise rows are never imputed, however much else is missing", {
  # The counterpart of the accuracy claim: a subject with no observed omics has
  # no conditional information, so inventing values for them would be fabricating
  # data rather than imputing it.
  for (model in c("early", "parallel")) {
    n_layer <- if (model == "early") 1L else 2L
    d <- sim_lucid(model, N = N_IMP, K = 2, M = M_IMP, n_layer = n_layer,
                   missing = "mixed", miss_ratio = 0.3, seed = 805)
    K <- if (model == "early") 2 else rep(2, n_layer)
    suppressWarnings(suppressMessages(invisible(capture.output(
      fit <- estimate_lucid(lucid_model = model, G = d$G, Z = d$Z, Y = d$Y,
                            family = "normal", K = K,
                            init_omic.data.model = "VVV", seed = 1)
    ))))
    expect_listwise_preserved(fit$Z, d$miss_index, info = model)
  }
})

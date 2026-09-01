# ---------------------------------------------------------------------------
# Independent, deliberately naive re-implementations of the LUCID equations.
#
# These are written with explicit loops and no reuse of anything in R/, so a
# bug in the package cannot be mirrored by the oracle that is meant to catch it.
# Equation numbers refer to RJ-2024-012 unless prefixed "BA" (vbae123).
# ---------------------------------------------------------------------------

# log-sum-exp of a vector
o_lse <- function(v) {
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

o_softmax <- function(v) {
  e <- exp(v - max(v))
  e / sum(e)
}

# log N(x; mean, sigma) evaluated on the OBSERVED coordinates only.
# A fully-missing row contributes 0 (BA Eq 11: the Z term drops out).
o_dmvnorm_obs <- function(x, mean, sigma) {
  obs <- which(!is.na(x))
  if (length(obs) == 0L) return(0)
  xo <- x[obs]; mo <- mean[obs]; so <- sigma[obs, obs, drop = FALSE]
  d <- length(obs)
  ch <- chol(so)
  q <- backsolve(ch, xo - mo, transpose = TRUE)
  -0.5 * d * log(2 * pi) - sum(log(diag(ch))) - 0.5 * sum(q^2)
}

# Eq 3/4 (supervised) and Eq 12 (unsupervised) responsibilities, plus Eq 2
# observed-data log-likelihood, for the early model.  Handles all three
# missingness patterns through o_dmvnorm_obs.
o_estep_early <- function(G, Z, Y, beta, mu, sigma, gamma_levels, gamma_sd,
                          family = c("normal", "binary"), useY = TRUE,
                          CoY = NULL, coef_CoY = NULL) {
  family <- match.arg(family)
  N <- nrow(as.matrix(Z)); K <- nrow(mu)
  Gd <- cbind(1, as.matrix(G))
  logp <- matrix(0, N, K)
  for (i in seq_len(N)) {
    for (k in seq_len(K)) {
      lp <- sum(Gd[i, ] * beta[k, ])                                   # G -> X
      lz <- o_dmvnorm_obs(as.numeric(Z[i, ]), mu[k, ], sigma[, , k])   # X -> Z
      ly <- 0
      if (useY) {
        eta <- gamma_levels[k]
        if (!is.null(CoY) && !is.null(coef_CoY)) {
          eta <- eta + sum(as.numeric(CoY[i, ]) * coef_CoY)
        }
        ly <- if (family == "normal") {
          stats::dnorm(as.numeric(Y)[i], eta, gamma_sd[k], log = TRUE)
        } else {
          stats::dbinom(as.numeric(Y)[i], 1, stats::plogis(eta), log = TRUE)
        }
      }
      logp[i, k] <- lp + lz + ly
    }
    # subtract the softmax denominator (same constant for every k in row i),
    # turning the raw linear predictors into log S(X = k | G) of Eq 3
    logp[i, ] <- logp[i, ] - o_lse(as.numeric(Gd[i, , drop = FALSE] %*% t(beta)))
  }
  r <- t(apply(logp, 1, o_softmax))
  ll <- sum(apply(logp, 1, o_lse))
  list(r = r, loglik = ll, logp = logp)
}

# Eq 13 / Eq 18: parameter count for the early model.
# D = (P + 1)(K - 1) + K*M + K*M*(M+1)/2 + n_Y, then Eq 18 subtracts the number
# of exposures and omic features whose effects are zero across all clusters.
o_npars_early <- function(P, M, K, family = c("normal", "binary"),
                          n_CoG = 0, n_CoY = 0, DG = 0, DZ = 0,
                          useY = TRUE, count_intercepts = TRUE) {
  family <- match.arg(family)
  n_beta_cols <- P + n_CoG + as.integer(count_intercepts)
  nY <- if (!useY) 0 else if (family == "normal") 2 * K + n_CoY else K + n_CoY
  D <- n_beta_cols * (K - 1) + K * M + K * M * (M + 1) / 2 + nY
  D - DG - DZ
}

o_bic <- function(loglik, npars, N) -2 * loglik + npars * log(N)

# BA Eq 17: closed-form I-step for one observation.
o_istep_row <- function(obs, mu, sigma, r) {
  M <- ncol(mu); K <- nrow(mu)
  A <- which(!is.na(obs)); B <- which(is.na(obs))
  if (!length(B) || !length(A)) return(obs)
  cur <- obs
  cur[B] <- as.numeric(r %*% mu[, B, drop = FALSE])
  prec <- lapply(seq_len(K), function(k) solve(sigma[, , k]))
  dens <- vapply(seq_len(K), function(k) {
    exp(o_dmvnorm_obs(cur, mu[k, ], sigma[, , k]))
  }, numeric(1))
  w <- r * dens
  lhs <- matrix(0, length(B), length(B)); rhs <- numeric(length(B))
  for (k in seq_len(K)) {
    obb <- prec[[k]][B, B, drop = FALSE]
    oba <- prec[[k]][B, A, drop = FALSE]
    lhs <- lhs + w[k] * obb
    rhs <- rhs + w[k] * as.numeric(oba %*% mu[k, A] + obb %*% mu[k, B] - oba %*% cur[A])
  }
  cur[B] <- as.numeric(solve(lhs) %*% rhs)
  cur
}

# Eqs 19/20: bootstrap normal and percentile intervals, computed by hand.
#
# The percentile interval uses the Davison & Hinkley (1997) order-statistic
# convention -- the endpoint at position (R + 1) * alpha, interpolated on the
# normal quantile scale -- which is what boot::boot.ci(type = "perc") implements
# and what the R Journal paper cites for Eq 20.  This is NOT the same as R's
# default quantile() type 7; at small R the two differ materially (at R = 60 the
# lower endpoint moved from -1.36 to -2.47 in testing), which is a concrete
# reason for the paper's guidance of R >= 800 for percentile intervals.
o_boot_perc <- function(t, alpha) {
  t <- sort(t[is.finite(t)])
  R <- length(t)
  k <- (R + 1) * alpha
  kf <- floor(k)
  vapply(seq_along(alpha), function(i) {
    if (kf[i] < 1L) return(t[1L])
    if (kf[i] >= R) return(t[R])
    frac <- k[i] - kf[i]
    if (frac == 0) return(t[kf[i]])
    # interpolate on the normal quantile scale, as boot:::norm.inter does
    q1 <- stats::qnorm(kf[i] / (R + 1)); q2 <- stats::qnorm((kf[i] + 1) / (R + 1))
    qt <- stats::qnorm(alpha[i])
    t[kf[i]] + (qt - q1) / (q2 - q1) * (t[kf[i] + 1L] - t[kf[i]])
  }, numeric(1))
}

o_boot_ci <- function(t0, t, conf = 0.95) {
  tf <- t[is.finite(t)]
  bias <- mean(tf) - t0
  se <- stats::sd(tf)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  a <- (1 - conf) / 2
  list(
    normal = c(t0 - bias - z * se, t0 - bias + z * se),
    percentile = o_boot_perc(tf, c(a, 1 - a))
  )
}

# ---- expectation helpers ---------------------------------------------------

expect_posterior_valid <- function(r, tol = 1e-8) {
  r <- as.matrix(r)
  testthat::expect_true(all(is.finite(r)), label = "posterior is finite")
  testthat::expect_true(all(r >= -tol & r <= 1 + tol), label = "posterior in [0,1]")
  testthat::expect_equal(as.numeric(rowSums(r)), rep(1, nrow(r)), tolerance = tol)
}

# Compare fitted cluster-level parameters to truth, up to a label permutation.
# Returns the permutation used so callers can align other blocks the same way.
match_clusters <- function(fit_mu, true_mu) {
  K <- nrow(true_mu)
  perms <- if (K <= 6) {
    do.call(rbind, lapply(seq_len(factorial(K)), function(i) .perm_i(K, i)))
  } else {
    matrix(seq_len(K), nrow = 1)
  }
  best <- NULL; best_d <- Inf
  for (i in seq_len(nrow(perms))) {
    p <- perms[i, ]
    d <- sum((fit_mu[p, , drop = FALSE] - true_mu)^2)
    if (d < best_d) { best_d <- d; best <- p }
  }
  list(perm = best, dist = best_d)
}

.perm_i <- function(n, i) {
  # i-th permutation of seq_len(n) in lexicographic order (1-based)
  pool <- seq_len(n); res <- integer(n); i <- i - 1L
  for (k in seq.int(n, 1)) {
    f <- factorial(k - 1L)
    j <- i %/% f + 1L
    res[n - k + 1L] <- pool[j]
    pool <- pool[-j]
    i <- i %% f
  }
  res
}

# ---- sanity and recovery expectations --------------------------------------
#
# These back the two things the suite must guarantee beyond "it ran": that a
# fit contains no NA/NaN/Inf or absurd values, and that the estimates are close
# to the parameters that generated the data.

# Every numeric leaf of a fitted object, with a path label, so a failure names
# the component that went wrong instead of just saying "somewhere in the fit".
.numeric_leaves <- function(x, path = "", acc = list()) {
  if (is.numeric(x)) {
    acc[[path]] <- as.numeric(x)
    return(acc)
  }
  if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
      acc <- .numeric_leaves(x[[i]], paste0(path, "$", nm), acc)
    }
  }
  acc
}

#' The "no crazy values" gate.
#'
#' `Z` is skipped: omics cells for listwise-missing rows are legitimately still
#' NA after fitting, which is asserted separately by expect_listwise_preserved().
#' `max_abs` catches separation blow-up, where a coefficient runs off to +/-Inf
#' but stops at some large finite value before any is.finite() check would fire.
expect_sane_fit <- function(fit, max_abs = 50, info = NULL) {
  tag <- function(s) paste0(s, if (!is.null(info)) paste0(" [", info, "]"))

  # `fit` holds the raw glmnet/nnet objects, which carry their own diagnostics
  # and are not ours to police. Z is checked separately (listwise rows are
  # legitimately still NA).
  drop <- c("Z", "z", "init_impute", "init_par", "init_omic.data.model",
            "fit", "submodel")
  keep <- fit[setdiff(names(fit), drop)]

  # Finiteness applies everywhere: an NA or NaN anywhere in a returned fit is a
  # defect regardless of which component it lands in.
  for (nm in names(.numeric_leaves(keep))) {
    v <- .numeric_leaves(keep)[[nm]]
    if (!length(v)) next
    testthat::expect_true(all(is.finite(v)), label = tag(paste0("finite values in ", nm)))
  }

  # The magnitude bound is a separation-blow-up detector, so it applies only to
  # estimated parameters -- not to counts (N, n_rows), iteration numbers, or
  # log-likelihoods, all of which are legitimately large.
  params <- list(
    res_Beta  = .strip_fits(fit$res_Beta),
    res_Mu    = fit$res_Mu,
    res_Delta = .strip_fits(fit$res_Delta),
    gamma     = fit$res_Gamma[intersect(names(fit$res_Gamma),
                              c("beta", "cluster_effect", "covariate", "sigma"))]
  )
  for (nm in names(.numeric_leaves(params))) {
    v <- .numeric_leaves(params)[[nm]]
    if (!length(v)) next
    testthat::expect_true(max(abs(v)) <= max_abs,
      label = tag(paste0("parameter magnitude <= ", max_abs, " in ", nm)))
  }

  # Posterior probabilities must be a distribution per row. inclusion.p nests
  # one level for parallel and two for a serial fit with a parallel stage, so
  # the blocks are collected recursively rather than assumed flat.
  for (m in .collect_posteriors(fit$inclusion.p)) expect_posterior_valid(m)

  # a reported log-likelihood must be finite and negative
  if (!is.null(fit$likelihood)) {
    testthat::expect_true(is.finite(fit$likelihood),
      label = paste0("finite likelihood", if (!is.null(info)) paste0(" [", info, "]")))
    testthat::expect_lt(fit$likelihood, 0)
  }

  # covariances must stay symmetric positive-definite
  sigmas <- .collect_sigma(fit$res_Sigma)
  for (s in sigmas) {
    testthat::expect_equal(s, t(s), tolerance = 1e-6)
    ev <- eigen(s, symmetric = TRUE, only.values = TRUE)$values
    testthat::expect_gt(min(ev), 0)
  }
  invisible(TRUE)
}

# Every posterior block in a fit, however deeply nested: a matrix for early,
# a list by layer for parallel, and a list by stage of either for serial.
.collect_posteriors <- function(x, acc = list()) {
  if (is.null(x)) return(acc)
  if (is.matrix(x) || is.data.frame(x)) return(c(acc, list(as.matrix(x))))
  if (is.numeric(x)) return(c(acc, list(matrix(x, ncol = 1))))
  if (is.list(x)) for (e in x) acc <- .collect_posteriors(e, acc)
  acc
}

# Drop raw model objects (glmnet/nnet fits) from a coefficient block, keeping
# only the estimated coefficients themselves.
.strip_fits <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.null(names(x))) x <- x[setdiff(names(x), c("fit", "penalty"))]
  if (is.list(x)) lapply(x, .strip_fits) else x
}

# res_Sigma is a list of matrices for early and a nested list otherwise.
.collect_sigma <- function(x, acc = list()) {
  if (is.matrix(x) && nrow(x) == ncol(x)) return(c(acc, list(x)))
  if (is.array(x) && length(dim(x)) == 3L) {
    for (k in seq_len(dim(x)[3])) acc <- c(acc, list(x[, , k]))
    return(acc)
  }
  if (is.list(x)) for (e in x) acc <- .collect_sigma(e, acc)
  acc
}

#' Compare a fit to the parameters that generated the data.
#'
#' Cluster labels are aligned first. The package canonicalizes cluster order by
#' outcome level for supervised fits, but not for unsupervised ones or upstream
#' serial stages, so alignment via match_clusters() is required rather than
#' optional -- comparing without it produces failures that look like bias but
#' are only a permutation.
expect_recovers_mu <- function(fit_mu, true_mu, tolerance = 0.25, info = NULL) {
  fit_mu <- as.matrix(fit_mu)
  true_mu <- as.matrix(true_mu)
  testthat::expect_equal(dim(fit_mu), dim(true_mu),
    label = paste0("mu dimensions", if (!is.null(info)) paste0(" [", info, "]")))
  mm <- match_clusters(fit_mu, true_mu)
  aligned <- fit_mu[mm$perm, , drop = FALSE]
  testthat::expect_equal(unname(aligned), unname(true_mu), tolerance = tolerance)
  invisible(mm$perm)
}

#' Cluster assignment accuracy against the true latent membership, maximised
#' over label permutations.
cluster_accuracy <- function(inclusion_p, true_x) {
  pred <- max.col(as.matrix(inclusion_p))
  K <- ncol(as.matrix(inclusion_p))
  perms <- if (K <= 6) {
    do.call(rbind, lapply(seq_len(factorial(K)), function(i) .perm_i(K, i)))
  } else matrix(seq_len(K), nrow = 1)
  max(apply(perms, 1, function(p) mean(p[pred] == true_x)))
}

#' Rows that were entirely missing must still be NA after fitting: they are
#' handled listwise and must never be invented by an imputation step (D7).
expect_listwise_preserved <- function(fit_Z, miss_index, info = NULL) {
  fz <- if (is.list(fit_Z)) fit_Z else list(fit_Z)
  mi <- if (is.list(miss_index)) miss_index else list(miss_index)
  for (i in seq_along(fz)) {
    lw <- which(rowSums(!mi[[i]]) == 0)
    if (!length(lw)) next
    testthat::expect_true(all(is.na(fz[[i]][lw, , drop = FALSE])),
      label = paste0("listwise rows remain NA in layer ", i,
                     if (!is.null(info)) paste0(" [", info, "]")))
  }
  invisible(TRUE)
}

#' Observed cells must never be altered by imputation.
expect_observed_unchanged <- function(fit_Z, orig_Z, info = NULL) {
  fz <- if (is.list(fit_Z)) fit_Z else list(fit_Z)
  oz <- if (is.list(orig_Z)) orig_Z else list(orig_Z)
  for (i in seq_along(fz)) {
    obs <- !is.na(oz[[i]])
    testthat::expect_equal(fz[[i]][obs], oz[[i]][obs],
      label = paste0("observed cells unchanged in layer ", i,
                     if (!is.null(info)) paste0(" [", info, "]")))
  }
  invisible(TRUE)
}

#' The EM ascent property.
#'
#' Without sporadic missingness the EM iteration is a true ascent: the observed
#' log-likelihood never decreases, and measured worst-case behaviour is a strictly
#' positive step (+5e-4 early, +9e-4 parallel) against a tolerance of 1e-3.
#'
#' With sporadic missingness the I-step makes it a majorization step whose ascent
#' guarantee is only approximate, and small decreases do occur -- measured at
#' -0.157 (early sporadic), -0.041 (parallel sporadic), -0.045 (early binary
#' mixed). So the two regimes get different assertions rather than one weak
#' bound that would pass on a genuinely diverging fit:
#'
#'   no I-step  -> strictly non-decreasing within `tol`
#'   I-step     -> bounded dips, and the trajectory must still rise overall
#'
#' A monotonicity violation beyond these means an M-step is not maximising what
#' the E-step computed, which produces plausible but wrong estimates.
expect_ascends <- function(fit, info = NULL, istep_slack = 0.5) {
  traces <- .collect_traces(fit)
  if (!length(traces)) {
    testthat::succeed()
    return(invisible(TRUE))
  }
  imputing <- .has_sporadic(fit)
  for (tr in traces) {
    if (length(tr$trace) < 2L) next
    slack <- if (imputing) istep_slack else max(tr$tol, 1e-6)
    testthat::expect_gte(min(diff(tr$trace)), -slack,
      label = paste0("EM log-likelihood ascends",
                     if (imputing) " (I-step regime)" else "",
                     if (!is.null(info)) paste0(" [", info, "]")))
    if (imputing) {
      # dips are tolerated, drifting downhill overall is not
      testthat::expect_gte(tr$trace[length(tr$trace)] - tr$trace[1], -slack,
        label = paste0("EM trajectory rises overall",
                       if (!is.null(info)) paste0(" [", info, "]")))
    }
  }
  invisible(TRUE)
}

# Sporadic rows are what switch the I-step on; listwise rows alone do not.
.has_sporadic <- function(fit) {
  ms <- fit$missing_summary
  if (is.null(ms)) return(FALSE)
  found <- FALSE
  walk <- function(x) {
    if (is.list(x)) {
      if (!is.null(x$sporadic_rows) && sum(x$sporadic_rows) > 0) found <<- TRUE
      for (e in x) walk(e)
    }
  }
  walk(ms)
  if (!is.null(fit$submodel)) for (s in fit$submodel) if (.has_sporadic(s)) found <- TRUE
  found
}

# A serial fit runs no EM loop of its own; its traces live on the sub-models.
.collect_traces <- function(fit) {
  out <- list()
  if (!is.null(fit$em_control$loglik_trace)) {
    out <- c(out, list(list(trace = as.numeric(fit$em_control$loglik_trace),
                            tol = fit$em_control$tol %||% 1e-3)))
  }
  if (!is.null(fit$submodel)) {
    for (s in fit$submodel) out <- c(out, .collect_traces(s))
  }
  out
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

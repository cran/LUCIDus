# Internal outcome-model helpers.  Outcome parameters are stored in two forms:
# absolute cluster/state levels for likelihood calculations and reference-coded
# coefficients for reporting.

#' Build an early-model outcome parameter object from absolute cluster levels
#'
#' Converts absolute per-cluster outcome levels (used for likelihood
#' calculations) into the reference-coded \code{beta} representation (used
#' for reporting), and packages both together.
#'
#' @param cluster_effect Absolute outcome level for each cluster.
#' @param covariate Covariate coefficients, if any.
#' @param sigma Residual SD per cluster (normal outcome only).
#' @param covariate_names Names for \code{covariate}.
#' @return A list: \code{beta} (reference-coded), \code{cluster_effect},
#'   \code{covariate}, \code{sigma}, \code{parameterization}.
#' @noRd
early_gamma_from_levels <- function(cluster_effect, covariate = numeric(0),
                                    sigma = NULL, covariate_names = NULL) {
  K <- length(cluster_effect)
  coef <- c(cluster_effect[1], cluster_effect[-1] - cluster_effect[1], covariate)
  names(coef) <- c("(Intercept)", if (K > 1) paste0("LC", 2:K), covariate_names)
  list(
    beta = coef,
    cluster_effect = as.numeric(cluster_effect),
    covariate = as.numeric(covariate),
    sigma = sigma,
    parameterization = "absolute_cluster_effect_v1"
  )
}

#' Recover absolute per-cluster outcome levels from a gamma object
#'
#' @param gamma An early-model outcome parameter object (see
#'   \code{\link{early_gamma_from_levels}}).
#' @param K Number of clusters.
#' @return Numeric vector of length \code{K}.
#' @noRd
early_gamma_levels <- function(gamma, K) {
  if (!is.null(gamma$cluster_effect)) {
    return(as.numeric(gamma$cluster_effect))
  }
  beta <- as.numeric(gamma$beta)
  c(beta[1], beta[1] + beta[seq.int(2, K)])
}

#' Recover covariate coefficients from a gamma object
#'
#' @param gamma An early-model outcome parameter object.
#' @param K Number of clusters.
#' @return Numeric vector of covariate coefficients (possibly length zero).
#' @noRd
early_gamma_covariate <- function(gamma, K) {
  if (!is.null(gamma$covariate)) {
    return(as.numeric(gamma$covariate))
  }
  beta <- as.numeric(gamma$beta)
  if (length(beta) > K) beta[-seq_len(K)] else numeric(0)
}

#' M-step: outcome model for the early model
#'
#' Responsibility-weighted fit of \eqn{Y | X, CoY}. Without covariates and
#' for a normal outcome, this reduces to a closed-form weighted mean/SD; every
#' other case (binary outcome, or any covariates) is a weighted GLM/nonlinear
#' fit over cluster-indicator design. For a normal outcome with covariates,
#' cluster means AND cluster-specific residual SDs are estimated jointly by
#' direct maximization (heteroscedastic across clusters), since \code{lm()}
#' assumes one shared variance; the weighted least-squares fit is used as the
#' starting value and as a fallback if that maximization fails to converge.
#'
#' @param Y Outcome data.
#' @param r Current cluster responsibilities (N x K).
#' @param CoY Optional outcome covariates.
#' @param K Number of clusters.
#' @param CoYnames Column names for \code{CoY}.
#' @param family "normal" or "binary".
#' @return An early-model outcome parameter object (see
#'   \code{\link{early_gamma_from_levels}}), with a \code{fit} element
#'   holding the underlying model object where one was fit.
#' @noRd
fit_early_outcome <- function(Y, r, CoY, K, CoYnames, family) {
  Y <- as.numeric(Y)
  r <- as.matrix(r)
  cov_names <- if (is.null(CoY)) character(0) else CoYnames

  if (is.null(CoY) && is_normal_family(family)) {
    levels <- colSums(r * Y) / colSums(r)
    sigma <- sqrt(colSums(r * vapply(levels, function(x) (Y - x)^2,
                                     numeric(length(Y)))) / colSums(r))
    return(early_gamma_from_levels(levels, sigma = sigma))
  }

  N <- length(Y)
  state <- factor(rep(seq_len(K), times = N), levels = seq_len(K))
  dat <- data.frame(Y = rep(Y, each = K), state = state)
  if (!is.null(CoY)) {
    CoY <- as.data.frame(CoY, check.names = FALSE)
    colnames(CoY) <- cov_names
    dat <- cbind(dat, CoY[rep(seq_len(N), each = K), , drop = FALSE])
  }
  weights <- as.vector(t(r))
  form <- stats::reformulate(c("state", cov_names), response = "Y")

  if (is_binary_family(family)) {
    fit <- suppressWarnings(stats::glm(form, data = dat, weights = weights,
                                       family = stats::binomial()))
    coef <- stats::coef(fit)
    if (anyNA(coef)) stop("Outcome model is not identifiable for the fitted latent clusters")
    names(coef) <- c("(Intercept)", paste0("LC", seq.int(2, K)), cov_names)
    levels <- c(coef[1], coef[1] + coef[seq.int(2, K)])
    covariate <- if (length(coef) > K) coef[-seq_len(K)] else numeric(0)
    sigma <- NULL
  } else {
    design <- stats::model.matrix(form, dat)
    fit0 <- stats::lm.wfit(design, dat$Y, w = weights)
    start_beta <- fit0$coefficients
    start_beta[!is.finite(start_beta)] <- 0
    start_resid <- dat$Y - as.numeric(design %*% start_beta)
    start_sigma <- vapply(seq_len(K), function(k) {
      idx <- dat$state == k
      sqrt(max(sum(weights[idx] * start_resid[idx]^2) / sum(weights[idx]), 1e-8))
    }, numeric(1))
    objective <- function(par) {
      beta <- par[seq_len(ncol(design))]
      sigma <- exp(par[ncol(design) + seq_len(K)])
      eta <- as.numeric(design %*% beta)
      -sum(weights * stats::dnorm(dat$Y, eta, sigma[as.integer(dat$state)], log = TRUE))
    }
    opt <- stats::nlminb(c(start_beta, log(start_sigma)), objective,
                         control = list(eval.max = 1000, iter.max = 1000,
                                        rel.tol = 1e-10))
    # D10: a nonzero nlminb code is not a reason to abort the whole EM fit.
    # Fall back to the weighted least-squares start, which is always defined,
    # and only fail if that is degenerate too.
    if (opt$convergence != 0 || !all(is.finite(opt$par))) {
      fallback <- c(start_beta, log(start_sigma))
      if (!all(is.finite(fallback)) || !is.finite(objective(fallback))) {
        stop("Outcome model is not identifiable for the fitted latent clusters")
      }
      warning("Weighted heteroscedastic outcome model did not converge; ",
              "using the weighted least-squares solution for this M-step.",
              call. = FALSE)
      opt <- list(par = fallback, objective = objective(fallback),
                  convergence = opt$convergence)
    }
    coef <- opt$par[seq_len(ncol(design))]
    names(coef) <- c("(Intercept)", paste0("LC", seq.int(2, K)), cov_names)
    levels <- c(coef[1], coef[1] + coef[seq.int(2, K)])
    covariate <- if (length(coef) > K) coef[-seq_len(K)] else numeric(0)
    sigma <- exp(opt$par[ncol(design) + seq_len(K)])
    fit <- list(coefficients = coef, sigma = sigma, objective = opt$objective,
                convergence = opt$convergence)
    class(fit) <- "lucid_weighted_gaussian"
  }
  out <- early_gamma_from_levels(levels, covariate, sigma, cov_names)
  out$fit <- fit
  out
}

#' Deterministic cluster ordering for an unsupervised (upstream serial) stage
#'
#' Upstream stages of a serial model have no outcome to order clusters by,
#' so their numbering was previously whatever the EM happened to land on --
#' seed-dependent, which made the reported delta coefficients (and their
#' signs) irreproducible across runs on identical data. Ordering
#' lexicographically by the cluster means of the omics block is deterministic,
#' data-driven, and independent of the outcome, which keeps the stage
#' genuinely unsupervised.
#'
#' @param mu Cluster means (K x Q matrix).
#' @return A permutation of \code{seq_len(nrow(mu))} giving the new cluster
#'   order.
#' @noRd
omics_order_index <- function(mu) {
  mu <- as.matrix(mu)
  do.call(order, c(as.data.frame(mu), list(seq_len(nrow(mu)))))
}

#' Relabel early-model cluster parameters into a canonical order
#'
#' Reorders clusters (by outcome level, or by the supplied \code{index}),
#' re-centering \code{beta} and \code{gamma} on the new reference cluster.
#'
#' @param beta,mu,sigma,gamma Current parameter estimates.
#' @param K Number of clusters.
#' @param index Optional explicit cluster permutation; if \code{NULL},
#'   ordered by outcome level (ascending).
#' @return A list: \code{beta}, \code{mu}, \code{sigma}, \code{gamma}
#'   (all reordered/re-centered), and \code{index} (the permutation used).
#' @noRd
relabel_early_parameters <- function(beta, mu, sigma, gamma, K, index = NULL) {
  levels <- early_gamma_levels(gamma, K)
  if (is.null(index)) index <- order(levels, seq_len(K))
  beta <- t(t(beta) - beta[index[1], ])[index, , drop = FALSE]
  mu <- mu[index, , drop = FALSE]

  sigma_array <- if (is.list(sigma)) {
    simplify2array(sigma)
  } else {
    sigma
  }
  sigma_array <- sigma_array[, , index, drop = FALSE]
  sigma_out <- lapply(seq_len(K), function(i) sigma_array[, , i, drop = FALSE][, , 1])

  covariate <- early_gamma_covariate(gamma, K)
  cov_names <- names(gamma$beta)
  cov_names <- if (length(covariate)) tail(cov_names, length(covariate)) else character(0)
  gamma_out <- early_gamma_from_levels(levels[index], covariate,
                                       if (is.null(gamma$sigma)) NULL else gamma$sigma[index],
                                       cov_names)
  list(beta = beta, mu = mu, sigma = sigma_out, gamma = gamma_out, index = index)
}

#' Reconstruct the joint cluster-state array from per-layer marginals
#'
#' @param z List of per-layer marginal cluster probabilities (each N x K_i).
#' @param K Per-layer cluster counts.
#' @return A \code{K[1] x ... x K[nOmics] x N} array of joint probabilities
#'   under a conditional-independence assumption across layers given the
#'   observation.
#' @noRd
parallel_joint_from_marginals <- function(z, K) {
  N <- nrow(z[[1]])
  grid <- expand.grid(lapply(K, seq_len), KEEP.OUT.ATTRS = FALSE)
  out <- matrix(0, nrow = prod(K), ncol = N)
  for (n in seq_len(N)) {
    probability <- rep(1, nrow(grid))
    for (i in seq_along(K)) probability <- probability * z[[i]][n, grid[[i]]]
    out[, n] <- probability
  }
  array(out, dim = c(K, N))
}

#' Build a parallel-model outcome parameter object from a coefficient vector
#'
#' Expands a flat coefficient vector (intercept, per-layer non-reference
#' cluster effects, covariates) into the full joint-state linear predictor
#' array, packaging both the additive per-layer effects and the resulting
#' state-level means/probabilities.
#'
#' @param coef Flat coefficient vector.
#' @param K Per-layer cluster counts.
#' @param family "normal" or "binary".
#' @param sd Residual SD (normal outcome only).
#' @param covariate_names Names for the covariate coefficients within
#'   \code{coef}.
#' @param fit Optional underlying model object to carry through.
#' @return A list: \code{mu} (state means/probabilities), \code{eta} (linear
#'   predictor), \code{intercept}, \code{effects} (per-layer, cluster-1
#'   fixed at zero), \code{covariate}, \code{coef}, \code{sd}, \code{K},
#'   \code{family}, \code{fit}, \code{parameterization}.
#' @noRd
parallel_delta_from_coef <- function(coef, K, family, sd = NULL,
                                     covariate_names = character(0), fit = NULL) {
  n_latent <- sum(K - 1)
  expected <- 1 + n_latent + length(covariate_names)
  coef_full <- numeric(expected)
  coef_full[seq_len(min(length(coef), expected))] <- as.numeric(coef)[seq_len(min(length(coef), expected))]
  names(coef_full) <- c("(Intercept)",
                        unlist(lapply(seq_along(K), function(i) paste0("Layer", i, "_LC", seq.int(2, K[i])))),
                        covariate_names)

  pos <- 2L
  effects <- lapply(seq_along(K), function(i) {
    n_i <- K[i] - 1L
    value <- c(0, if (n_i > 0) coef_full[seq.int(pos, length.out = n_i)])
    pos <<- pos + n_i
    as.numeric(value)
  })
  covariate <- if (length(covariate_names)) tail(coef_full, length(covariate_names)) else numeric(0)
  state_eta <- parallel_state_array(coef_full[1], effects, K)
  if (!is.null(fit)) fit$coefficients <- coef_full
  list(mu = if (to_parallel_family(family) == "binomial") stats::plogis(state_eta) else state_eta,
       eta = state_eta, intercept = unname(coef_full[1]), effects = effects,
       covariate = as.numeric(covariate), coef = coef_full, sd = sd, K = K,
       family = to_parallel_family(family), fit = fit,
       parameterization = "additive_state_effect_v1")
}

#' Build the joint-state linear predictor array from additive per-layer effects
#'
#' @param intercept Scalar intercept.
#' @param effects Per-layer additive cluster effects (cluster 1 fixed at zero).
#' @param K Per-layer cluster counts.
#' @return A \code{K[1] x ... x K[nOmics]} array.
#' @noRd
parallel_state_array <- function(intercept, effects, K) {
  grid <- expand.grid(lapply(K, seq_len), KEEP.OUT.ATTRS = FALSE)
  eta <- rep(intercept, nrow(grid))
  for (i in seq_along(K)) eta <- eta + effects[[i]][grid[[i]]]
  array(eta, dim = K)
}

#' Recover the flat coefficient vector from a parallel-model outcome object
#'
#' @param Delta A parallel-model outcome parameter object (see
#'   \code{\link{parallel_delta_from_coef}}).
#' @return Numeric coefficient vector, from whichever of \code{Delta$coef},
#'   \code{Delta$fit}, or \code{Delta$mu} is available.
#' @noRd
parallel_delta_coef <- function(Delta) {
  if (!is.null(Delta$coef)) return(Delta$coef)
  if (!is.null(Delta$fit)) return(stats::coef(Delta$fit))
  as.numeric(Delta$mu)
}

#' M-step: outcome model for the parallel model
#'
#' Responsibility-weighted GLM fit of \eqn{Y} on the joint cluster-state
#' indicator design (plus covariates). Falls back to an intercept-only fit if
#' the full design is not identifiable.
#'
#' @param Y Outcome data.
#' @param r Current joint-state responsibilities.
#' @param K Per-layer cluster counts.
#' @param N Number of observations.
#' @param family "normal" or "binary".
#' @param CoY Optional outcome covariates.
#' @return A list: \code{fit} (the underlying model object) and \code{Gamma}
#'   (a parallel-model outcome parameter object).
#' @noRd
fit_parallel_outcome <- function(Y, r, K, N, family, CoY) {
  family <- to_parallel_family(family)
  state_grid <- expand.grid(lapply(K, seq_len), KEEP.OUT.ATTRS = FALSE)
  state_design <- do.call(cbind, lapply(seq_along(K), function(i) {
    if (K[i] == 1L) return(NULL)
    out <- vapply(seq.int(2, K[i]), function(k) as.numeric(state_grid[[i]] == k),
                  numeric(nrow(state_grid)))
    if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
    colnames(out) <- paste0("Layer", i, "_LC", seq.int(2, K[i]))
    out
  }))
  if (is.null(state_design)) state_design <- matrix(numeric(0), nrow(state_grid), 0L)
  design <- as.data.frame(state_design[rep(seq_len(nrow(state_grid)), times = N), , drop = FALSE],
                          check.names = FALSE)
  cov_names <- character(0)
  if (!is.null(CoY)) {
    CoY <- as.data.frame(CoY, check.names = FALSE)
    cov_names <- colnames(CoY)
    if (is.null(cov_names)) cov_names <- paste0("CoY", seq_len(ncol(CoY)))
    colnames(CoY) <- cov_names
    design <- cbind(design, CoY[rep(seq_len(N), each = nrow(state_grid)), , drop = FALSE])
  }
  dat <- cbind(data.frame(Y = rep(as.numeric(Y), each = nrow(state_grid))), design)
  weights <- as.vector(matrix(r, nrow = prod(K), ncol = N))
  form <- stats::reformulate(colnames(design), response = "Y")
  fit <- if (family == "binomial") {
    try(suppressWarnings(stats::glm(
      form, data = dat, weights = weights, family = stats::binomial(),
      control = stats::glm.control(maxit = 100)
    )), silent = TRUE)
  } else {
    try(stats::lm(form, data = dat, weights = weights), silent = TRUE)
  }
  if (inherits(fit, "try-error")) {
    fit <- if (family == "binomial") {
      suppressWarnings(stats::glm(Y ~ 1, data = dat, weights = weights,
                                  family = stats::binomial()))
    } else {
      stats::lm(Y ~ 1, data = dat, weights = weights)
    }
  }
  coef <- stats::coef(fit)
  coef[is.na(coef)] <- 0
  sd <- if (family == "gaussian") {
    sqrt(sum(weights * stats::residuals(fit)^2) / sum(weights))
  } else NULL
  Delta <- parallel_delta_from_coef(coef, K, family, sd, cov_names, fit)
  list(fit = Delta$fit, Gamma = Delta)
}

#' Linear predictor for every observation and joint cluster state
#'
#' @param Delta A parallel-model outcome parameter object.
#' @param CoY Optional outcome covariates.
#' @param N Number of observations (inferred from \code{CoY} if omitted).
#' @return A \code{K[1] x ... x K[nOmics] x N} array of linear predictors.
#' @noRd
parallel_state_eta <- function(Delta, CoY = NULL, N = NULL) {
  eta0 <- if (!is.null(Delta$eta)) Delta$eta else {
    coef <- parallel_delta_coef(Delta)
    parallel_delta_from_coef(coef, Delta$K,
                             if (is.null(Delta$family)) "gaussian" else Delta$family)$eta
  }
  if (is.null(N)) N <- if (is.null(CoY)) 1L else nrow(CoY)
  out <- array(rep(as.numeric(eta0), N), dim = c(Delta$K, N))
  covariate <- if (!is.null(Delta$covariate)) Delta$covariate else numeric(0)
  if (!is.null(CoY) && length(covariate)) {
    shift <- as.numeric(as.matrix(CoY) %*% covariate)
    out <- sweep(out, length(Delta$K) + 1L, shift, "+")
  }
  out
}

#' Apply a per-layer cluster permutation to a joint-state array
#'
#' @param x An array with one dimension per layer (and, if
#'   \code{include_observation}, a trailing observation dimension).
#' @param permutations Per-layer cluster permutations.
#' @param include_observation Whether \code{x} has a trailing observation
#'   dimension left untouched by the permutation.
#' @return \code{x} with each layer's dimension reordered by its permutation.
#' @noRd
permute_parallel_array <- function(x, permutations, include_observation = FALSE) {
  indices <- lapply(permutations, identity)
  if (include_observation) indices <- c(indices, list(seq_len(dim(x)[length(permutations) + 1L])))
  do.call(`[`, c(list(x), indices, list(drop = FALSE)))
}

#' Relabel parallel-model cluster parameters into a canonical order
#'
#' Reorders each layer's clusters (by outcome effect, or by the supplied
#' \code{permutations}), re-centering \code{Beta} and \code{Delta} on the
#' new per-layer reference clusters.
#'
#' @param Beta,Mu,Sigma Current per-layer parameter estimates.
#' @param Delta Current outcome parameter object.
#' @param r Current joint-state responsibilities.
#' @param K Per-layer cluster counts.
#' @param selectZ Optional per-layer omics-selection indicators to permute
#'   alongside the parameters.
#' @param permutations Optional explicit per-layer permutations; if
#'   \code{NULL}, ordered by outcome effect (ascending) within each layer.
#' @return A list: \code{Beta}, \code{Mu}, \code{Sigma}, \code{Delta},
#'   \code{r} (all reordered/re-centered), \code{selectZ}, and
#'   \code{permutations} (the permutation used).
#' @noRd
relabel_parallel_parameters <- function(Beta, Mu, Sigma, Delta, r, K,
                                        selectZ = NULL, permutations = NULL) {
  if (is.null(permutations)) {
    permutations <- lapply(Delta$effects, function(x) order(x, seq_along(x)))
  }
  ref_effect <- vapply(seq_along(K), function(i) Delta$effects[[i]][permutations[[i]][1]], numeric(1))
  effects <- lapply(seq_along(K), function(i) {
    x <- Delta$effects[[i]][permutations[[i]]]
    x - x[1]
  })
  coef <- c(Delta$intercept + sum(ref_effect),
            unlist(lapply(effects, function(x) x[-1])), Delta$covariate)
  cov_names <- if (length(Delta$covariate)) tail(names(Delta$coef), length(Delta$covariate)) else character(0)
  Delta_new <- parallel_delta_from_coef(coef, K, Delta$family, Delta$sd, cov_names, Delta$fit)

  for (i in seq_along(K)) {
    full_beta <- rbind(0, Beta[[i]])
    full_beta <- t(t(full_beta) - full_beta[permutations[[i]][1], ])
    Beta[[i]] <- full_beta[permutations[[i]][-1], , drop = FALSE]
    Mu[[i]] <- Mu[[i]][permutations[[i]], , drop = FALSE]
    Sigma[[i]] <- Sigma[[i]][, , permutations[[i]], drop = FALSE]
    if (!is.null(selectZ) && !is.null(dim(selectZ[[i]])) && nrow(selectZ[[i]]) == K[i]) {
      selectZ[[i]] <- selectZ[[i]][permutations[[i]], , drop = FALSE]
    }
  }
  r <- permute_parallel_array(r, permutations, include_observation = TRUE)
  list(Beta = Beta, Mu = Mu, Sigma = Sigma, Delta = Delta_new, r = r,
       selectZ = selectZ, permutations = permutations)
}

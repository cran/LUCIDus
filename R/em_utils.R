# EM-machinery utilities shared by the fitting entry points: parameter
# initializers, log-likelihood/penalty math, and small array helpers.

#' Marginal cluster-inclusion probability for one omics layer
#'
#' Collapses the joint responsibility array down to one layer's marginal
#' probabilities, row by row.
#'
#' @param r Joint responsibility array (one dimension per layer, plus a
#'   trailing observation dimension).
#' @param N Number of observations.
#' @param layer Which layer's margin to extract.
#' @return An \code{N x K[layer]} matrix of marginal probabilities.
#' @noRd
compute_res_r <- function(r, N, layer) {
  r_margin <- t(sapply(1:N, function(j) {
    marginSums(lastInd(r,j), margin = layer)
  }))
  return(r_margin)
}


#' Observed-data log-likelihood from a joint log-likelihood array
#'
#' Sums \code{safe_log_sum_exp()} over the latent-state dimension(s) for
#' each observation.
#'
#' @param log_joint A joint log-likelihood array: a matrix (early model,
#'   observations as rows) or a multi-dimensional array (parallel model,
#'   observations as the last dimension).
#' @param observations "rows" (early) or "last" (parallel): which dimension
#'   indexes observations.
#' @return The scalar total observed-data log-likelihood.
#' @noRd
observed_loglik <- function(log_joint, observations = c("rows", "last")) {
  observations <- match.arg(observations)
  d <- dim(log_joint)
  if (is.null(d) || length(d) < 2L) {
    stop("log_joint must contain latent-state and observation dimensions")
  }

  if (observations == "rows") {
    if (length(d) != 2L) {
      stop("row-oriented log_joint must be a matrix")
    }
    return(sum(apply(log_joint, 1L, safe_log_sum_exp)))
  }

  n <- d[length(d)]
  sum(vapply(seq_len(n), function(i) {
    safe_log_sum_exp(as.numeric(lastInd(log_joint, i)))
  }, numeric(1)))
}

#' Observed-data log-likelihood for the parallel model
#'
#' Thin wrapper around \code{\link{observed_loglik}} for the parallel
#' model's observation-last array layout.
#'
#' @param Estep_array Joint log-likelihood array from \code{\link{Estep}}.
#' @param Estep_r Unused; accepted for a consistent call signature.
#' @return The scalar total observed-data log-likelihood.
#' @noRd
cal_loglik <- function(Estep_array, Estep_r = NULL) {
  observed_loglik(Estep_array, observations = "last")
}

#' Total penalty value for the early model's penalized log-likelihood
#'
#' Sums the L1 penalty on exposure coefficients (excluding covariates), the
#' L1 penalty on omics means, and the graphical-lasso penalty on omics
#' precision matrices -- so \code{penalized loglik = loglik - this}.
#'
#' @param beta,mu,sigma Current parameter estimates.
#' @param rho_g,rho_mu,rho_cov \code{Rho_G}, \code{Rho_Z_Mu}, \code{Rho_Z_Cov}.
#' @param dim_g_exposure Number of true exposure columns in \code{beta}
#'   (remaining columns are covariates, never penalized).
#' @return The scalar total penalty.
#' @noRd
early_penalty_value <- function(beta, mu, sigma, rho_g, rho_mu, rho_cov,
                                dim_g_exposure) {
  beta_penalty <- if (rho_g == 0) 0 else {
    rho_g * sum(abs(beta[, 1L + seq_len(dim_g_exposure), drop = FALSE]))
  }
  mu_penalty <- if (rho_mu == 0) 0 else rho_mu * sum(abs(mu))
  precision_penalty <- if (rho_cov == 0) 0 else {
    rho_cov * sum(vapply(seq_len(dim(sigma)[3]), function(k) {
      sum(abs(safe_solve(sigma[, , k])))
    }, numeric(1)))
  }
  beta_penalty + mu_penalty + precision_penalty
}

#' Log-sum-exp, as used by the parallel model's exposure E-step
#'
#' Equivalent to \code{safe_log_sum_exp()} but without its edge-case
#' handling for all-infinite input; kept as the separate implementation
#' \code{\link{f_GtoX}} already calls.
#'
#' @param vec A numeric vector.
#' @return The scalar log-sum-exp of \code{vec}.
#' @noRd
LogSumExp <- function(vec) {
  max_vec <- max(vec)
  trick <- max_vec + log(sum(exp(vec - max_vec)))
  return(trick)
}



#' Extract the n-th slice along an array's last dimension
#'
#' @param x An array.
#' @param n Index along the last dimension.
#' @return The n-th slice, as an array with one fewer dimension than
#'   \code{x}.
#' @noRd
lastInd <- function(x, n){
  d <- dim(x)
  d.new <- d[-length(d)]
  block.size <- prod(d.new)
  res <- x[(block.size * (n - 1) + 1):(block.size * n)]
  array(res, dim = d.new)
}



#' Validate a K (number-of-clusters) specification
#'
#' Every element of \code{K} must be a whole number >= 2 -- a single-cluster
#' model has no latent structure to estimate. Called from every fitting and
#' tuning entry point (\code{lucid()}, \code{est_lucid()}, \code{estimate_lucid()},
#' \code{tune_lucid_auxi()}) so a bad \code{K} is always reported at the call
#' the user actually made, not from deep inside the EM machinery.
#'
#' @param K Integer, or vector/list of integers for a "parallel"/"serial" model.
#' @return Invisibly \code{NULL}; called for its \code{stop()} side effect.
#' @noRd
check_K <- function(K) {
  for(x in K) {
    if(x != as.integer(x)) {
      stop("K should be a vector of integer")
    }
  }
  if(min(K) < 2) {
    stop("each element in K should be greater or equal than 2")
  }
}




#' Random initial cluster means, per layer (parallel model)
#'
#' @param K Per-layer cluster counts.
#' @param nZ Per-layer number of omics features.
#' @return A list, one \code{K[i] x nZ[i]} matrix per layer.
#' @noRd
initialize_Mu <- function(K, nZ) {
  nOmics <- length(K)
  Mu <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # For parallel model: K clusters x M features
    Mu[[i]] <- matrix(runif(K[i] * nZ[i], min = -1, max = 1),
                      nrow = K[i],  # K clusters in rows
                      ncol = nZ[i]) # M features in columns
  }
  return(Mu)
}



#' Identity-matrix initial cluster covariances, per layer (parallel model)
#'
#' @param K Per-layer cluster counts.
#' @param nZ Per-layer number of omics features.
#' @return A list, one \code{nZ[i] x nZ[i] x K[i]} array of identity
#'   matrices per layer.
#' @noRd
initialize_Sigma <- function(K, nZ) {
  nOmics <- length(K)
  Sigma <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # each element is an array of nZ[i] x nZ[i] x K[i]
    Sigma[[i]] <- array(0, dim = c(nZ[i], nZ[i], K[i]))
    for(j in 1:K[i]) {
      Sigma[[i]][, , j] <- diag(nZ[i])
    }
  }
  return(Sigma)
}





#' Model-based initial cluster means/covariances, per layer (parallel model)
#'
#' Runs \code{mclust::Mclust()} per layer (on rows with an observed value in
#' that layer) to get data-driven starting values, rather than the
#' random/identity initializers.
#'
#' @param K Per-layer cluster counts.
#' @param Z Omics data, one matrix per layer.
#' @param modelNames Per-layer \code{mclust} geometric model, or \code{NA}
#'   to let \code{mclust} choose.
#' @param na_pattern Per-layer missingness info from \code{check_na()}.
#' @return A list: \code{Mu}, \code{Sigma} (per layer), \code{z} (initial
#'   responsibilities, expanded back to all \code{N} rows), and
#'   \code{modelNames} (the geometric model actually used per layer, needed
#'   downstream since the M-step cannot re-derive an auto-selected choice).
#' @noRd
initialize_Mu_Sigma <- function(K, Z, modelNames, na_pattern) {
  nOmics <- length(K)
  Mu <- vector(mode = "list", length = nOmics)
  Sigma <- vector(mode = "list", length = nOmics)
  z <- vector(mode = "list", length = nOmics)
  # The geometric model each layer ended up with. When `modelNames` is NULL the
  # caller asked mclust to choose, and the choice has to be reported back: the
  # M-step calls mclust::mstep(), which needs a concrete model name and cannot
  # re-derive it. Mirrors `model.best <- mclust.fit$modelName` on the early path.
  resolved <- character(nOmics)
  for(i in 1:nOmics) {
    # Exclude all-missing rows for layer-specific GMM fit.
    idx_obs <- na_pattern[[i]]$indicator_na != 3
    temp_fit <- Mclust(data = Z[[i]][idx_obs, ],
                       G = K[i],
                       modelNames = modelNames[i])
    if (is.null(temp_fit)) {
      stop("mclust failed for omics layer ", i,
           " - set init_omic.data.model to `NULL` to conduct automatic model selection")
    }
    resolved[i] <- temp_fit$modelName
    # Transpose mean matrix for parallel model - means should be K x M
    Mu[[i]] <- t(temp_fit$parameters$mean)  # Transpose to get K x M
    Sigma[[i]] <- temp_fit$parameters$variance$sigma
    # Expand initial responsibilities back to full N rows so downstream
    # multinomial regressions align with G.
    z_full <- matrix(1 / K[i], nrow = nrow(Z[[i]]), ncol = K[i])
    z_full[idx_obs, ] <- temp_fit$z
    z[[i]] <- z_full
  }
  return(list(Mu = Mu,
              Sigma = Sigma,
              z = z,
              modelNames = resolved))
}



#' Initial outcome parameter object from initial responsibilities (parallel model)
#'
#' @param K Per-layer cluster counts.
#' @param CoY Optional outcome covariates.
#' @param family "normal" or "binary".
#' @param z Per-layer initial responsibilities (from
#'   \code{\link{initialize_Mu_Sigma}}).
#' @param Y Outcome data.
#' @return A parallel-model outcome parameter object (see
#'   \code{\link{parallel_delta_from_coef}}).
#' @noRd
initialize_Delta <- function(K, CoY, family = c("normal", "binary"),
                             z, Y) {
  joint_r <- parallel_joint_from_marginals(z, K)
  fit_parallel_outcome(
    Y = Y, r = joint_r, K = K, N = nrow(z[[1]]),
    family = to_parallel_family(family), CoY = CoY
  )$Gamma
}




#' Expand a flat coefficient vector into a joint-state array
#'
#' @param K Per-layer cluster counts.
#' @param mu Flat coefficient vector: intercept followed by each layer's
#'   non-reference cluster effects.
#' @return A \code{K[1] x ... x K[nOmics]} array (see
#'   \code{\link{parallel_state_array}}).
#' @noRd
vec_to_array <- function(K, mu) {
  pos <- 2L
  effects <- lapply(K, function(k) {
    n <- k - 1L
    value <- c(0, if (n > 0L) mu[seq.int(pos, length.out = n)])
    pos <<- pos + n
    value
  })
  parallel_state_array(mu[1], effects, K)
}

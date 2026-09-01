# Parallel-integration E-step machinery. Parameters throughout this file:
# K (number of latent clusters per layer), Beta (G -> X coefficients),
# Mu/Sigma (X -> Z coefficients), Delta (X -> Y coefficients).

#' Log posterior cluster probability for one omics layer given exposures
#'
#' @param G Exposure design matrix.
#' @param Beta_matrix Non-reference-cluster exposure coefficients for this
#'   layer (cluster 1's row of zeros is prepended internally).
#' @return An \code{N x K} matrix of log \eqn{p(X_j = k \mid G)}.
#' @noRd
f_GtoX <- function(G, Beta_matrix) {
  N <- nrow(G)
  Beta_matrix <- rbind(rep(0, ncol(Beta_matrix)),
                       Beta_matrix)
  xb <- cbind(rep(1, N), G) %*% t(Beta_matrix)
  xb_LSE <- apply(xb, 1, LogSumExp)
  return(xb - xb_LSE)
}


#' Log likelihood for one omics layer given its latent cluster
#'
#' Evaluates the Gaussian mixture density on observed coordinates only
#' (\code{\link{log_dmvnorm_observed}}) -- previously any \code{NA} in
#' \code{Z} reached \code{mclust::dmvnorm} directly and raised
#' \code{"NA/NaN/Inf in foreign function call"}, which made prediction on
#' sporadically-missing omics impossible.
#'
#' @param Z Omics data for this layer.
#' @param Mu_matrix Cluster means (K x M).
#' @param Sigma_matrix Cluster covariances (M x M x K).
#' @return An \code{N x K} matrix of log \eqn{p(Z_j \mid X_j = k)}.
#' @noRd
f_XtoZ <- function(Z, Mu_matrix, Sigma_matrix) {
  N <- nrow(Z)
  M <- ncol(Z)
  K <- nrow(Mu_matrix)  # Number of clusters
  
  # Validate dimensions
  if (ncol(Mu_matrix) != M) {
    stop(sprintf("Dimension mismatch: Z has %d features but Mu_matrix has %d features", 
                M, ncol(Mu_matrix)))
  }
  if (!identical(dim(Sigma_matrix)[1:2], c(M, M))) {
    stop(sprintf("Invalid Sigma dimensions: Expected first two dimensions to be c(%d,%d), got c(%s)", 
                M, M, paste(dim(Sigma_matrix)[1:2], collapse=",")))
  }
  if (dim(Sigma_matrix)[3] != K) {
    stop(sprintf("Invalid Sigma dimensions: Expected %d clusters, got %d", 
                K, dim(Sigma_matrix)[3]))
  }
  
  XtoZ <- matrix(rep(0, N * K), nrow = N)
  for (i in 1:K) {
    mean_i <- Mu_matrix[i, ]
    sigma_i <- Sigma_matrix[, , i]
    
    # Additional validation
    if (length(mean_i) != M) {
      stop(sprintf("Mean vector for cluster %d has wrong length: expected %d, got %d", 
                  i, M, length(mean_i)))
    }
    
    # Ensure sigma is symmetric and positive definite
    sigma_i <- (sigma_i + t(sigma_i))/2  # Ensure symmetry
    eigen_values <- eigen(sigma_i, symmetric = TRUE)$values
    if (any(eigen_values <= 0)) {
      # Add small constant to diagonal if not positive definite
      sigma_i <- sigma_i + diag(1e-6, M)
    }

    # Evaluate on observed coordinates only.  Previously any NA in Z reached
    # mclust::dmvnorm and raised "NA/NaN/Inf in foreign function call", which
    # made prediction on sporadically-missing omics impossible.
    tmp_ll <- try({
      log_dmvnorm_observed(Z, mean = mean_i, sigma = sigma_i)
    }, silent = TRUE)

    # Check for errors
    if (inherits(tmp_ll, "try-error")) {
      stop(sprintf("Error in cluster %d: %s", i, attr(tmp_ll, "condition")$message))
    }
    XtoZ[, i] <- tmp_ll
  }
  return(XtoZ)
}


#' Log likelihood of the outcome given every joint cluster state
#'
#' @param Y Outcome data.
#' @param Delta Outcome parameter object (see
#'   \code{\link{parallel_delta_from_coef}}).
#' @param family "normal" or "binary".
#' @param CoY Optional outcome covariates.
#' @return A \code{K[1] x ... x K[nOmics] x N} array of log
#'   \eqn{p(Y \mid X = \text{state})}.
#' @noRd
f_XtoY <- function(Y, Delta, family, CoY = NULL) {
  family <- to_parallel_family(family)
  if (!is.matrix(Y)) stop("Y should be a matrix")

  N <- nrow(Y)
  K <- Delta$K
  eta <- parallel_state_eta(Delta, CoY = CoY, N = N)
  y_long <- rep(as.numeric(Y), each = prod(K))

  if (family == "gaussian") {
    return(array(
      stats::dnorm(y_long, mean = as.numeric(eta), sd = Delta$sd, log = TRUE),
      dim = c(K, N)
    ))
  }

  array(
    stats::dbinom(y_long, size = 1, prob = as.numeric(stats::plogis(eta)),
                  log = TRUE),
    dim = c(K, N)
  )
}


#' E-step: joint log-likelihood over every cluster state (parallel model)
#'
#' Combines the exposure (\eqn{G \to X}), omics (\eqn{X \to Z}, one term
#' per layer), and, if \code{useY}, outcome (\eqn{X \to Y}) log-likelihoods
#' for every joint combination of per-layer cluster states.
#'
#' @param G,Z,Y Exposure, omics (list, one per layer), and outcome data.
#' @param Beta,Mu,Sigma Current per-layer parameter estimates.
#' @param Delta Current outcome parameter object.
#' @param family "normal" or "binary".
#' @param useY Whether to include the outcome term.
#' @param na_pattern Per-layer missingness info from \code{check_na()}; a
#'   row with no observed value in a layer contributes 0 for that layer.
#' @param CoY Optional outcome covariates.
#' @return A \code{K[1] x ... x K[nOmics] x N} array of joint log-likelihoods.
#' @noRd
Estep <- function(G, Z, Y, Beta, Mu, Sigma, Delta, family, useY, na_pattern,
                  CoY = NULL) {
  family <- to_parallel_family(family)
  N <- nrow(G)
  K <- Delta$K
  res <- array(rep(0, prod(K) * N),
               dim = c(K, N))
  
  # Build layer-wise Z|X log-likelihood matrices on the full N rows.
  # Rows with all-missing omics values contribute 0 (no Z information).
  f2_all <- lapply(seq_along(K), function(i) {
    layer_res <- matrix(0, nrow = N, ncol = K[i])
    idx_obs <- na_pattern[[i]]$indicator_na != 3
    if (any(idx_obs)) {
      layer_res[idx_obs, ] <- f_XtoZ(
        Z = Z[[i]][idx_obs, , drop = FALSE],
        Mu_matrix = Mu[[i]],
        Sigma_matrix = Sigma[[i]]
      )
    }
    layer_res
  })


  f1 <- lapply(seq_along(K), function(i) f_GtoX(G, Beta[[i]]))
  states <- expand.grid(lapply(K, seq_len), KEEP.OUT.ATTRS = FALSE)
  res_matrix <- matrix(0, nrow = prod(K), ncol = N)
  for (s in seq_len(nrow(states))) {
    for (a in seq_along(K)) {
      state <- states[s, a]
      res_matrix[s, ] <- res_matrix[s, ] + f1[[a]][, state] + f2_all[[a]][, state]
    }
  }
  if (useY) {
    res_matrix <- res_matrix + matrix(
      f_XtoY(Y = Y, Delta = Delta, family = family, CoY = CoY),
      nrow = prod(K), ncol = N
    )
  }
  array(res_matrix, dim = c(K, N))
}



#' Normalize the E-step joint log-likelihood array to responsibilities
#'
#' @param Estep_array Joint log-likelihood array from \code{\link{Estep}}.
#' @param K Per-layer cluster counts.
#' @param N Number of observations.
#' @return A \code{K[1] x ... x K[nOmics] x N} array of joint
#'   responsibilities, each observation's slice summing to one.
#' @noRd
Estep_to_r <- function(Estep_array, K, N) {

  out <- matrix(0, nrow = prod(K), ncol = N)
  raw <- matrix(Estep_array, nrow = prod(K), ncol = N)
  for (i in seq_len(N)) out[, i] <- lse_vec(raw[, i])
  array(out, dim = c(K, N))
}

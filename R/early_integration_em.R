#' Randomly initialize a covariance matrix per cluster
#'
#' Draws each cluster's starting covariance from a uniform distribution over
#' a symmetric matrix, then stabilizes it.
#'
#' @param dimZ Number of omics features.
#' @param K Number of latent clusters.
#' @return A \code{dimZ x dimZ x K} array, one stabilized covariance matrix
#'   per cluster.
#' @noRd
gen_cov_matrices <- function(dimZ, K) {
  res <- array(rep(0, dimZ^2 * K), dim = c(dimZ, dimZ, K))
  for(i in 1:K) {
    x <- matrix(runif(dimZ^2, min = -0.5, max = 0.5), nrow = dimZ)
    x_sym <- t(x) %*% x
    # Ensure numerical stability of initial covariance
    res[,,i] <- check_and_stabilize_sigma(x_sym)
  }
  return(res)
}



#' E-step: per-observation, per-cluster joint log-likelihood (early model)
#'
#' Computes \eqn{\log p(G, Z, Y, X = k)} for each observation and cluster
#' \code{k} -- the exposure term (\eqn{G \to X}), the omics term
#' (\eqn{X \to Z}, evaluated on observed coordinates only, which reproduces
#' the list-wise likelihood partition of vbae123 Eq 11 for a fully-missing
#' row), and, if \code{useY}, the outcome term (\eqn{X \to Y}).
#'
#' @param beta,mu,sigma,gamma Current parameter estimates (exposure, omics
#'   mean/covariance, outcome).
#' @param G,Z,Y Exposure, omics, and outcome data.
#' @param family.list Outcome-family dispatch object from
#'   \code{\link{normal}}/\code{\link{binary}}.
#' @param K,N Number of clusters and observations.
#' @param useY Whether to include the outcome term.
#' @param ind.na Per-row missingness indicator from \code{check_na()}.
#' @param ... Passed to \code{family.list$f.pYgX}.
#' @return An \code{N x K} matrix of joint log-likelihoods.
#' @noRd
Estep_early <- function(beta,
                  mu,
                  sigma,
                  gamma = NULL,
                  G,
                  Z,
                  Y = NULL,
                  family.list,
                  K,
                  N,
                  useY,
                  ind.na, ...) {
  # initialize vectors for storing likelihood
  pXgG <- pZgX <- pYgX <- matrix(0, nrow = N, ncol = K)

  # Safe log-likelihood for G -> X
  xb <- cbind(rep(1, N), G) %*% t(beta)
  xb_lse <- apply(xb, 1, safe_log_sum_exp)
  pXgG <- sweep(xb, 1, xb_lse, "-")

  # Safe log-likelihood for X -> Z.
  # Rows are evaluated on their OBSERVED coordinates only, which is the correct
  # marginal Gaussian under the model.  Rows with no observed omics value
  # (ind.na == 3) contribute 0, reproducing the list-wise likelihood partition
  # of vbae123 Eq 11.  This makes training and prediction share one code path:
  # previously a row with any NA fell through to a fallback that returned NA,
  # collapsing the posterior to uniform.
  keep <- ind.na != 3
  for (i in 1:K) {
    # Stabilize covariance matrix
    sigma[,,i] <- check_and_stabilize_sigma(sigma[,,i])
    if (any(keep)) {
      pZgX[keep, i] <- log_dmvnorm_observed(
        Z[keep, , drop = FALSE],
        mean = mu[i, ],
        sigma = sigma[, , i]
      )
    }
  }

  # log-likelihood for X->Y
  if(useY){
    pYgX <- family.list$f.pYgX(Y, gamma, K = K, N = N, ...)
  }

  vec <- pXgG + pZgX + pYgX
  
  # Handle numerical issues
  vec[is.na(vec)] <- -Inf
  vec[is.infinite(vec) & vec > 0] <- .Machine$double.max / 2
  vec[is.infinite(vec) & vec < 0] <- -.Machine$double.max / 2
  
  return (vec)
}





#' M-step: exposure-to-cluster association (early model)
#'
#' Fits the multinomial logistic \eqn{G \to X} model, penalized (LASSO via
#' \code{glmnet}) if \code{selectG}, unpenalized (\code{nnet::multinom})
#' otherwise. Falls back to the unpenalized fit if the LASSO fit fails.
#'
#' @param G Exposure (and covariate) design matrix.
#' @param r Current cluster responsibilities (N x K).
#' @param selectG Whether to apply exposure selection.
#' @param penalty LASSO penalty (\code{Rho_G}).
#' @param dimG,dimCoG Number of exposure and covariate columns.
#' @param K Number of clusters.
#' @return A \code{K x (dimG + dimCoG + 1)} coefficient matrix, cluster 1
#'   as reference.
#' @noRd
Mstep_G <- function(G, r, selectG, penalty, dimG, dimCoG, K) {
  new.beta <- matrix(0, nrow = K, ncol = (dimG + dimCoG + 1))
  if(selectG){
    if(dimG < 2) {
      stop("At least 2 exposure variables are needed for variable selection")
    }
    penalty.factor <- c(rep(1, dimG), rep(0, dimCoG))
    tryLasso <- try(glmnet::glmnet(
      x = as.matrix(G),
      y = as.matrix(r),
      family = "multinomial",
      lambda = penalty,
      penalty.factor = penalty.factor
    ), silent = TRUE)
    if(inherits(tryLasso, "try-error")) {
      warning("Lasso failed, using unpenalized estimation")
      beta.multilogit <- nnet::multinom(as.matrix(r) ~ as.matrix(G), trace = FALSE)
      new.beta[-1, ] <- coef(beta.multilogit)
    } else {
      new.beta[, 1] <- tryLasso$a0
      new.beta[, -1] <- t(matrix(unlist(lapply(tryLasso$beta, function(x) x[,1])), ncol = K))
    }
  }
  else{
    beta.multilogit <- nnet::multinom(as.matrix(r) ~ as.matrix(G), trace = FALSE)
    new.beta[-1, ] <- coef(beta.multilogit)
  }
  return(new.beta)
}



#' M-step: cluster-specific omics mean and covariance (early model)
#'
#' Penalized (via \code{\link{penalized_cluster_block}}) or unpenalized
#' (\code{mclust::mstep}) Gaussian mixture M-step, restricted to rows with at
#' least one observed omics value.
#'
#' @param Z Omics data.
#' @param r Current cluster responsibilities.
#' @param selectZ Whether to apply omics selection.
#' @param penalty.mu,penalty.cov LASSO and graphical-lasso penalties
#'   (\code{Rho_Z_Mu}, \code{Rho_Z_Cov}).
#' @param model.name \code{mclust} geometric model for the covariance.
#' @param K Number of clusters.
#' @param ind.na Per-row missingness indicator; rows with no observed value
#'   (\code{== 3}) are excluded.
#' @param mu Current cluster means, used as the starting point for the
#'   penalized estimator.
#' @return A list with \code{mu} (\code{K x Q} matrix) and \code{sigma}
#'   (\code{Q x Q x K} array).
#' @noRd
Mstep_Z <- function(Z, r, selectZ, penalty.mu, penalty.cov,
                    model.name, K, ind.na, mu) {
  dz <- Z[ind.na != 3, ]
  dr <- r[ind.na != 3, ]
  Q <- ncol(Z)
  new_sigma <- array(0, dim = c(Q, Q, K))
  new_mu <- matrix(0, nrow = K, ncol = Q)
  if(selectZ) {
    for(k in 1:K) {
      blk <- penalized_cluster_block(
        z = dz, r = dr, k = k, mu_k = mu[k, ],
        penalty.mu = penalty.mu, penalty.cov = penalty.cov
      )
      new_mu[k, ] <- blk$mu
      new_sigma[, , k] <- blk$sigma
    }
  }
  else{
    z.fit <- mclust::mstep(modelName = model.name, data = dz, z = dr)
    new_mu <- t(z.fit$parameters$mean)
    new_sigma <- z.fit$parameters$variance$sigma
    
    # Ensure positive definiteness for all components
    for(k in 1:K) {
      new_sigma[, , k] <- check_and_stabilize_sigma(new_sigma[, , k])
    }
  }
  return(list(mu = new_mu, sigma = new_sigma))
}


#' Log-sum-exp, as a call-time wrapper around \code{safe_log_sum_exp}
#'
#' Defined as a wrapper rather than the alias \code{lse <- safe_log_sum_exp}:
#' an alias is evaluated when this file is sourced, which would force
#' \code{safe_log_sum_exp} to exist before this file in collation order. A
#' wrapper resolves the name at call time instead, so no load-order coupling
#' is needed.
#'
#' @param x A numeric vector or matrix; see \code{safe_log_sum_exp}.
#' @return See \code{safe_log_sum_exp}.
#' @noRd
lse <- function(x) safe_log_sum_exp(x)

#' Normalize a vector of log-likelihoods to a probability vector
#'
#' The log-sum-exp trick (\code{exp(vec - safe_log_sum_exp(vec))}), with
#' explicit fallbacks to the uniform distribution for \code{NA}, all-\code{-Inf},
#' or all-\code{Inf} input, and a final renormalization if floating-point
#' error leaves the result not quite summing to one.
#'
#' @param vec A vector of length \code{K} log-likelihoods for one observation.
#' @return A probability vector of length \code{K} summing to one.
#' @noRd
lse_vec <- function(vec) {
  # Check for invalid inputs
  if (any(is.na(vec))) {
    warning("NA values in probability calculation")
    return(rep(1/length(vec), length(vec)))
  }
  
  # Handle infinite values
  if (any(is.infinite(vec))) {
    if (all(is.infinite(vec) & vec < 0)) {
      warning("All negative infinite values in probability calculation")
      return(rep(1/length(vec), length(vec)))
    }
    if (all(is.infinite(vec) & vec > 0)) {
      warning("All positive infinite values in probability calculation")
      return(rep(1/length(vec), length(vec)))
    }
  }
  
  # Use safe normalization
  norm_vec <- safe_normalize(exp(vec - safe_log_sum_exp(vec)))
  
  # Final check for numerical stability
  if (any(is.na(norm_vec) | is.infinite(norm_vec))) {
    warning("Numerical instability in probability normalization")
    return(rep(1/length(vec), length(vec)))
  }
  
  # Ensure probabilities sum to 1 (within numerical precision)
  if (abs(sum(norm_vec) - 1) > sqrt(.Machine$double.eps)) {
    norm_vec <- norm_vec / sum(norm_vec)
  }
  
  return(norm_vec)
}

#' Closed-form L1-penalized mean update for one cluster (Zhou et al.)
#'
#' Soft-thresholding update for the omics mean of cluster \code{j}, driven
#' by the precision matrix \code{wi} (the graphical-lasso or empirical
#' precision), used by \code{\link{penalized_cluster_block}}.
#'
#' @param j Cluster index.
#' @param rho L1 penalty on the mean (\code{Rho_Z_Mu}).
#' @param z Observed omics rows (n x Q).
#' @param r Responsibilities for the same rows (n x K).
#' @param mu Current mean for cluster \code{j} (length Q).
#' @param wi Precision matrix for cluster \code{j} (Q x Q).
#' @return The updated mean for cluster \code{j}, length Q, with entries
#'   thresholded to exactly zero where the penalty dominates.
#' @noRd
est_mu <- function(j, rho, z, r, mu, wi){
  p <- ncol(z)
  res.mu <- numeric(p)
  for(x in 1:p) {
    q1 <- t(t(z) - mu) %*% wi[x, ]
    q2 <- q1 + wi[x, x] * z[, x] - wi[x, x] * (z[, x] - mu[x])
    sum_q2r <- sum(q2 * r[, j])
    
    if(abs(sum_q2r) <= rho) {
      res.mu[x] <- 0
    } else {
      a <- sum(r[, j] * rowSums(t(wi[x, ] * t(z))))
      b <- sum(r[, j]) * (sum(wi[x, ] * mu) - wi[x, x] * mu[x])
      t1 <- (a - b + rho) / (sum(r[, j]) * wi[x, x])
      t2 <- (a - b - rho) / (sum(r[, j]) * wi[x, x])
      
      res.mu[x] <- if(t1 < 0) t1 else t2
    }
  }
  return(res.mu)
}

#' Shared penalized M-step block for one cluster (Eqs 15-17)
#'
#' \code{mu} is updated by the L1-penalized estimator of Zhou et al.
#' (\code{\link{est_mu}}), driven by the graphical-lasso precision matrix;
#' \code{sigma} is the graphical-lasso covariance. Used by both the early
#' model (\code{Mstep_Z}) and the parallel model (\code{Mstep_XtoZ}) so that
#' penalized estimation is one estimator, not two. When \code{penalty.cov <= 0}
#' the empirical covariance is used directly rather than calling
#' \code{glasso::glasso()}, which would otherwise emit a rank warning on
#' every call -- hundreds of identical warnings per fit.
#'
#' @param z Observed omics rows (n x Q).
#' @param r Responsibilities for the same rows (n x K).
#' @param k Cluster index.
#' @param mu_k Current mean for cluster \code{k} (length Q).
#' @param penalty.mu,penalty.cov LASSO and graphical-lasso penalties.
#' @return A list with \code{mu} and \code{sigma} for cluster \code{k}.
#' @noRd
penalized_cluster_block <- function(z, r, k, mu_k, penalty.mu, penalty.cov) {
  z <- as.matrix(z)
  r <- as.matrix(r)
  Q <- ncol(z)
  w <- r[, k]
  w_sum <- sum(w)
  if (!is.finite(w_sum) || w_sum <= 0) {
    w <- rep(1, nrow(z))
    w_sum <- length(w)
  }

  # Empirical covariance about the current mean (Eq 17)
  Z_mu <- t(t(z) - mu_k)
  E_S <- crossprod(Z_mu * sqrt(w)) / w_sum

  # With penalty.cov == 0 the graphical lasso reduces to the empirical
  # covariance, but glasso() still emits a rank warning on every call -- which
  # fires once per cluster per layer per EM iteration (hundreds of identical
  # warnings).  Take the equivalent closed form directly instead.
  if (penalty.cov <= 0) {
    sigma_k <- E_S
    wi <- safe_solve(check_and_stabilize_sigma(E_S))
    mu_new <- est_mu(j = k, rho = penalty.mu, z = z, r = r,
                     mu = mu_k, wi = wi)
  } else {
    l_cov <- try(glasso::glasso(E_S, penalty.cov), silent = TRUE)
    if (inherits(l_cov, "try-error")) {
      warning("Glasso failed, using unpenalized estimation", call. = FALSE)
      sigma_k <- E_S
      mu_new <- colSums(w * z) / w_sum
    } else {
      sigma_k <- l_cov$w
      mu_new <- est_mu(j = k, rho = penalty.mu, z = z, r = r,
                       mu = mu_k, wi = l_cov$wi)
    }
  }

  list(mu = mu_new, sigma = check_and_stabilize_sigma(sigma_k))
}

####################Early Integration Estep####################

# Use uniform distributon to initialzie variance covariance matrices
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



# Calculate the log-likelihood of cluster assignment for each observation
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

  # Safe log-likelihood for X -> Z
  for (i in 1:K) {
    # Stabilize covariance matrix
    sigma[,,i] <- check_and_stabilize_sigma(sigma[,,i])
    if (any(ind.na != 3)) {
      z_ll <- try({
        mclust::dmvnorm(
          data = Z[ind.na != 3,],
          mean = mu[i,],
          sigma = sigma[,,i],
          log = TRUE
        )
      }, silent = TRUE)
      
      # Handle any errors in mvnorm calculation
      if (inherits(z_ll, "try-error")) {
        warning("Error in multivariate normal calculation, using fallback")
        # Use diagonal covariance as fallback
        diag_sigma <- diag(diag(sigma[,,i]))
        z_diag <- dnorm(
          Z[ind.na != 3,, drop = FALSE],
          mean = mu[i,],
          sd = sqrt(diag(diag_sigma)),
          log = TRUE
        )
        pZgX[ind.na != 3, i] <- rowSums(z_diag)
      } else {
        pZgX[ind.na != 3, i] <- z_ll
      }
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





# M-step to estimate the association between exposure and latent cluster
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
    ))
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



# M-step to estimate the parameters related to GMM for omics data
Mstep_Z <- function(Z, r, selectZ, penalty.mu, penalty.cov,
                    model.name, K, ind.na, mu) {
  dz <- Z[ind.na != 3, ]
  dr <- r[ind.na != 3, ]
  Q <- ncol(Z)
  new_sigma <- array(0, dim = c(Q, Q, K))
  new_mu <- matrix(0, nrow = K, ncol = Q)
  if(selectZ) {
    for(k in 1:K) {
      # Estimate E(S_k) for glasso
      Z_mu <- t(t(dz) - mu[k, ])
      E_S <- (matrix(colSums(dr[, k] * t(apply(Z_mu, 1, function(x) x %*% t(x)))), 
                     Q, Q)) / sum(dr[, k])
      
      # Use glasso with error handling
      l_cov <- try(glasso::glasso(E_S, penalty.cov))
      
      if(inherits(l_cov, "try-error")) {
        warning("Glasso failed, using unpenalized estimation")
        # Fallback to sample covariance
        new_sigma[, , k] <- E_S
        new_mu[k, ] <- colSums(dr[, k] * dz) / sum(dr[, k])
      } else {
        new_sigma[, , k] <- l_cov$w
        new_mu[k, ] <- est_mu(
          j = k,
          rho = penalty.mu,
          z = dz,
          r = dr,
          mu = mu[k, ],
          wi = l_cov$wi
        )
      }
      
      # Ensure positive definiteness
      new_sigma[, , k] <- check_and_stabilize_sigma(new_sigma[, , k])
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


# calculate the log-sum-exp
lse <- safe_log_sum_exp

# use the log-sum-exp trick to normalize a vector to probability
# @param vec  a vector of length K
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

# M-step: obtain sparse mean for omics data via LASSO penalty
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

#' E-step for early integration model
#'
#' @param G Exposure matrix
#' @param Z Omics data matrix
#' @param Y Outcome vector
#' @param CoG Covariates for G (optional)
#' @param CoY Covariates for Y (optional)
#' @param params Current parameter estimates
#' @param family Outcome distribution family
#' @return List containing posterior probabilities and log-likelihood
#' @noRd
early_estep <- function(G, Z, Y, CoG = NULL, CoY = NULL, params, family) {
  n_obs <- nrow(G)
  K <- nrow(params$mu)
  
  # Initialize log-probabilities matrix
  log_probs <- matrix(0, n_obs, K)
  
  # G->X contribution
  for(k in 2:K) {  # First cluster is reference
    log_probs[,k] <- G %*% params$beta[k-1,]
    if(!is.null(CoG)) {
      log_probs[,k] <- log_probs[,k] + CoG %*% params$beta_cov[k-1,]
    }
  }
  
  # X->Z contribution
  for(k in 1:K) {
    centered <- scale(Z, center = params$mu[k,], scale = FALSE)
    sigma_inv <- safe_solve(params$sigma[,,k])
    log_det_sigma <- determinant(params$sigma[,,k], logarithm = TRUE)$modulus
    
    log_probs[,k] <- log_probs[,k] - 
      0.5 * rowSums((centered %*% sigma_inv) * centered) -
      0.5 * log_det_sigma -
      0.5 * ncol(Z) * log(2 * pi)
  }
  
  # X->Y contribution (if using Y)
  if(!is.null(params$gamma) && !is.null(Y)) {
    if(family == "normal") {
      for(k in 1:K) {
        y_mean <- params$gamma$beta[k]
        if(!is.null(CoY)) {
          y_mean <- y_mean + CoY %*% params$gamma$beta_cov
        }
        log_probs[,k] <- log_probs[,k] -
          0.5 * ((Y - y_mean)^2 / params$gamma$sigma^2) -
          0.5 * log(2 * pi * params$gamma$sigma^2)
      }
    } else if(family == "binary") {
      for(k in 1:K) {
        logit <- params$gamma$beta[k]
        if(!is.null(CoY)) {
          logit <- logit + CoY %*% params$gamma$beta_cov
        }
        prob <- 1 / (1 + exp(-logit))
        log_probs[,k] <- log_probs[,k] +
          Y * safe_log(prob) + (1 - Y) * safe_log(1 - prob)
      }
    }
  }
  
  # Convert to probabilities using log-sum-exp trick
  log_norm <- apply(log_probs, 1, safe_log_sum_exp)
  posterior <- exp(sweep(log_probs, 1, log_norm))
  
  # Ensure numerical stability
  posterior <- pmax(pmin(posterior, 1), 0)
  posterior <- posterior / rowSums(posterior)
  
  # Calculate log-likelihood
  loglik <- sum(log_norm)
  
  return(list(
    posterior = posterior,
    loglik = loglik
  ))
}

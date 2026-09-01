#' M-step: exposure-to-cluster association, per layer (parallel model)
#'
#' Fits one multinomial logistic \eqn{G \to X_j} model per omics layer, on
#' that layer's marginal responsibilities, penalized (LASSO via
#' \code{glmnet}, cross-validated if \code{penalty} is \code{NULL}) or
#' unpenalized. \code{glmnet}'s multinomial fit returns one coefficient
#' block per class with no reference; the fitted coefficients are rebased
#' onto class 1 to match the unpenalized branch's \code{K - 1}-row,
#' reference-coded shape (the softmax is invariant to which class is used as
#' the reference).
#'
#' @param G Exposure (and covariate) design matrix.
#' @param r Current joint-state responsibilities.
#' @param selectG Whether to apply exposure selection.
#' @param penalty LASSO penalty (\code{Rho_G}); cross-validated if
#'   \code{NULL} and \code{selectG}.
#' @param K Per-layer cluster counts.
#' @param N Number of observations.
#' @param dimG_exposure Number of true exposure columns in \code{G}
#'   (remaining columns are covariates, never penalized); defaults to all of
#'   \code{ncol(G)}.
#' @return A list: \code{Beta} (per-layer coefficient matrices), and
#'   \code{selectG}/\code{selectG_layer} (selection indicators, if
#'   \code{selectG}).
#' @noRd
Mstep_GtoX <- function(G, r, selectG, penalty, K, N, dimG_exposure = NULL) {
  nOmics <- length(K)
  dimG_total <- ncol(G)
  if(is.null(dimG_exposure)) {
    dimG_exposure <- dimG_total
  }
  dimG_exposure <- as.integer(dimG_exposure)
  if(is.na(dimG_exposure) || dimG_exposure < 1 || dimG_exposure > dimG_total) {
    stop("dimG_exposure should be between 1 and ncol(G)")
  }
  dimCoG <- dimG_total - dimG_exposure
  penalty.factor <- c(rep(1, dimG_exposure), rep(0, dimCoG))
  # store multinomial logistic regression model with corresponding coefficients
  fit <- vector(mode = "list", length = nOmics)
  Beta <- vector(mode = "list", length = nOmics)
  selectG_results <- vector(mode = "list", length = nOmics)
  
  for(i in 1:nOmics) {
    r_margin <- t(sapply(1:N, function(j) {
      marginSums(lastInd(r,j), margin = i)
    }))
    
    if(selectG) {
      if(dimG_exposure < 2) {
        stop("At least 2 exposure variables are needed for feature selection")
      }
      
      # Try multiple lambda values if penalty is NULL
      if(is.null(penalty)) {
        cv_fit <- cv.glmnet(as.matrix(G), 
                           as.matrix(r_margin),
                           family = "multinomial",
                           alpha = 1,
                           penalty.factor = penalty.factor)
        penalty <- cv_fit$lambda.min
      }
      
      # Fit lasso with optimal/provided penalty
      tryLasso <- try({
        fit_lasso <- glmnet(as.matrix(G),
                           as.matrix(r_margin),
                           family = "multinomial",
                           lambda = penalty,
                           penalty.factor = penalty.factor)
        
        # D2: glmnet's multinomial fit returns one coefficient block per
        # CLASS (K blocks), using a symmetric parameterization with no
        # reference.  f_GtoX() prepends a zero reference row, so storing all K
        # blocks produced K + 1 states and the E-step then read only the first
        # K columns of a (K+1)-normalized matrix, silently discarding ~1/(K+1)
        # of the posterior mass.  Rebase onto class 1 and keep exactly K - 1
        # rows, matching the unpenalized nnet::multinom branch.  The softmax is
        # invariant to which class is used as the reference.
        beta_coef <- coef(fit_lasso)
        beta_full <- do.call(rbind, lapply(beta_coef, function(x) x[, 1]))
        Beta[[i]] <- sweep(beta_full, 2L, beta_full[1L, ], "-")[-1L, , drop = FALSE]
        
        # Determine selected features based on the coefficient range across
        # ALL K clusters, including the zero reference row.
        beta_exposure <- rbind(0, Beta[[i]])[, 1 + seq_len(dimG_exposure), drop = FALSE]
        coef_range <- apply(beta_exposure, 2, function(x) diff(range(x)))
        selectG_results[[i]] <- as.logical(abs(coef_range) > 0.001)
        
        fit[[i]] <- fit_lasso
        TRUE  # Return TRUE to indicate success
      }, silent = TRUE)
      
      if(inherits(tryLasso, "try-error")) {
        warning(sprintf("Lasso failed for omics layer %d, falling back to unpenalized regression", i))
        # Fall back to unpenalized regression
        temp_fit <- suppressWarnings(nnet::multinom(r_margin ~ G, trace = FALSE))
        fit[[i]] <- temp_fit
        Beta[[i]] <- coef(temp_fit)
        selectG_results[[i]] <- rep(TRUE, dimG_exposure)
      }
    } else {
      # Unpenalized regression
      temp_fit <- suppressWarnings(nnet::multinom(r_margin ~ G, trace = FALSE))
      fit[[i]] <- temp_fit
      Beta[[i]] <- coef(temp_fit)
      selectG_results[[i]] <- rep(TRUE, dimG_exposure)
    }
    
    # Add column names
    if(!is.null(colnames(G))) {
      colnames(Beta[[i]])[-1] <- colnames(G)
    }
  }
  
  return(list(
    fit = fit,
    Beta = Beta,
    selectG = selectG_results,
    penalty = penalty
  ))
}

#' M-step: cluster-specific omics mean and covariance, per layer (parallel model)
#'
#' Per-layer Gaussian mixture M-step on that layer's marginal
#' responsibilities, restricted to rows with an observed value in that
#' layer; penalized via \code{\link{penalized_cluster_block}} (the same
#' estimator the early model uses) or unpenalized via \code{mclust::mstep}.
#'
#' @param Z Omics data, one matrix per layer.
#' @param r Current joint-state responsibilities.
#' @param K Per-layer cluster counts.
#' @param modelNames Per-layer \code{mclust} geometric model.
#' @param N Number of observations.
#' @param na_pattern Per-layer missingness info from \code{check_na()}.
#' @param selectZ Whether to apply omics selection.
#' @param penalty.mu,penalty.cov LASSO and graphical-lasso penalties.
#' @return A list: \code{Mu}, \code{Sigma} (per layer), and \code{selectZ}
#'   (per-layer selection indicators, if \code{selectZ}).
#' @noRd
Mstep_XtoZ <- function(Z, r, K, modelNames, N, na_pattern, selectZ = FALSE, penalty.mu = 0, penalty.cov = 0,
                       mu = NULL) {
  nOmics <- length(K)
  # store GMM model with corresponding model
  fit <- vector(mode = "list", length = nOmics)
  Mu <- vector(mode = "list", length = nOmics)
  Sigma <- vector(mode = "list", length = nOmics)
  selectZ_results <- vector(mode = "list", length = nOmics)
  
  for(i in 1:nOmics) {
    r_margin <- t(sapply(1:N, function(j) {
      marginSums(lastInd(r,j), margin = i)
    }))
    r_margin <- round(r_margin, digits = 8)
    
    if(selectZ) {
      # Feature selection for Z using graphical lasso
      idx_obs <- na_pattern[[i]]$indicator_na != 3
      Z_i <- Z[[i]][idx_obs, , drop = FALSE]
      
      # Estimate means and covariances for each cluster.  The penalized
      # estimator is iterative: est_mu() and the empirical covariance are both
      # taken about the CURRENT mean, so seed from the previous iteration's Mu
      # (falling back to the weighted mean on the first pass) rather than 0.
      Mu[[i]] <- if (!is.null(mu) && !is.null(mu[[i]]) &&
                     identical(dim(as.matrix(mu[[i]])), c(as.integer(K[i]), ncol(Z_i)))) {
        as.matrix(mu[[i]])
      } else {
        t(vapply(seq_len(K[i]), function(k) {
          w <- r_margin[idx_obs, k]
          if (!is.finite(sum(w)) || sum(w) <= 0) w <- rep(1, nrow(Z_i))
          colSums(w * Z_i) / sum(w)
        }, numeric(ncol(Z_i))))
      }
      Sigma[[i]] <- array(0, dim = c(ncol(Z_i), ncol(Z_i), K[i]))
      selectZ_results[[i]] <- matrix(FALSE, nrow = K[i], ncol = ncol(Z_i))
      
      # D3: use the SAME penalized estimator as the early model
      # (L1-penalized mean via est_mu() driven by the glasso precision matrix,
      # glasso covariance) instead of an independent soft-threshold.
      r_obs <- r_margin[idx_obs, , drop = FALSE]
      for(k in 1:K[i]) {
        blk <- penalized_cluster_block(
          z = Z_i, r = r_obs, k = k, mu_k = Mu[[i]][k, ],
          penalty.mu = penalty.mu, penalty.cov = penalty.cov
        )
        Mu[[i]][k, ] <- blk$mu
        Sigma[[i]][, , k] <- blk$sigma

        # D3: selection is driven by the MEAN only.  The previous condition
        # OR-ed in colSums(abs(Sigma)) > 0.001, which is always TRUE because a
        # valid covariance has a positive diagonal -- so every feature was
        # reported as selected even when its mean was shrunk to exactly zero.
        selectZ_results[[i]][k, ] <- abs(Mu[[i]][k, ]) > 0.001
      }
      
      fit[[i]] <- list(
        parameters = list(
          mean = t(Mu[[i]]),
          variance = list(sigma = Sigma[[i]])
        )
      )
    } else {
      # Standard GMM without feature selection
      temp_fit <- mstep(data = Z[[i]][na_pattern[[i]]$indicator_na != 3, ],
                       G = K[i],
                       z = r_margin[na_pattern[[i]]$indicator_na != 3, ],
                       modelName = modelNames[i])
      fit[[i]] <- temp_fit
      Mu[[i]] <- t(temp_fit$parameters$mean)  # Transpose to get K x M format
      Sigma[[i]] <- temp_fit$parameters$variance$sigma
      selectZ_results[[i]] <- matrix(TRUE, nrow = K[i], ncol = ncol(Z[[i]]))
    }
  }
  
  # Update log-likelihood with regularization penalties
  loglik <- 0
  if(selectZ) {
    for(i in 1:nOmics) {
      if(penalty.mu > 0) {
        loglik <- loglik - penalty.mu * sum(abs(Mu[[i]]))
      }
      if(penalty.cov > 0) {
        for(k in 1:K[i]) {
          loglik <- loglik - penalty.cov * sum(abs(Sigma[[i]][,,k]))
        }
      }
    }
  }
  
  return(list(
    fit = fit,
    Mu = Mu,
    Sigma = Sigma,
    selectZ = selectZ_results,
    loglik = loglik
  ))
}

#' M-step: outcome model (parallel model)
#'
#' Thin wrapper around \code{\link{fit_parallel_outcome}}.
#'
#' @param Y Outcome data.
#' @param r Current joint-state responsibilities.
#' @param K Per-layer cluster counts.
#' @param N Number of observations.
#' @param family "normal" or "binary".
#' @param CoY Optional outcome covariates.
#' @return See \code{fit_parallel_outcome()}.
#' @noRd
Mstep_XtoY <- function(Y, r, K, N, family, CoY) {
  fit_parallel_outcome(
    Y = Y, r = r, K = K, N = N,
    family = to_parallel_family(family), CoY = CoY
  )
}

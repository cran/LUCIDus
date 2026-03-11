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
        
        # Extract coefficients and convert to matrix
        beta_coef <- coef(fit_lasso)
        Beta[[i]] <- do.call(rbind, lapply(beta_coef, function(x) x[,1]))
        
        # Determine selected features based on coefficient range
        beta_exposure <- Beta[[i]][, 1 + seq_len(dimG_exposure), drop = FALSE]
        coef_range <- apply(beta_exposure, 2, function(x) diff(range(x)))
        selectG_results[[i]] <- as.logical(abs(coef_range) > 0.001)
        
        fit[[i]] <- fit_lasso
        TRUE  # Return TRUE to indicate success
      })
      
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

Mstep_XtoZ <- function(Z, r, K, modelNames, N, na_pattern, selectZ = FALSE, penalty.mu = 0, penalty.cov = 0) {
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
      
      # Estimate means and covariances for each cluster
      Mu[[i]] <- matrix(0, nrow = K[i], ncol = ncol(Z_i))
      Sigma[[i]] <- array(0, dim = c(ncol(Z_i), ncol(Z_i), K[i]))
      selectZ_results[[i]] <- matrix(FALSE, nrow = K[i], ncol = ncol(Z_i))
      
      for(k in 1:K[i]) {
        # Weighted mean estimation with L1 penalty
        weights <- r_margin[idx_obs, k]
        w_sum <- sum(weights)
        if (!is.finite(w_sum) || w_sum <= 0) {
          weights <- rep(1, nrow(Z_i))
          w_sum <- nrow(Z_i)
        }
        weighted_mean <- colSums(weights * Z_i) / w_sum
        
        # L1-penalized mean estimation
        if(penalty.mu > 0) {
          for(j in 1:ncol(Z_i)) {
            soft_threshold <- sign(weighted_mean[j]) * max(0, abs(weighted_mean[j]) - penalty.mu)
            Mu[[i]][k,j] <- soft_threshold
          }
        } else {
          Mu[[i]][k,] <- weighted_mean
        }
        
        # Graphical lasso for covariance estimation
        Z_centered <- scale(Z_i, center = Mu[[i]][k,], scale = FALSE)
        S <- cov.wt(Z_centered, wt = weights)$cov
        if(penalty.cov > 0) {
          gl_fit <- try({
            glasso_result <- glasso(S, rho = penalty.cov)
            Sigma[[i]][,,k] <- glasso_result$w
            TRUE  # Return TRUE to indicate success
          })
          
          if(inherits(gl_fit, "try-error")) {
            warning(sprintf("Graphical lasso failed for layer %d, cluster %d. Using regularized covariance.", i, k))
            Sigma[[i]][,,k] <- S + diag(penalty.cov, ncol(S))
          }
        } else {
          Sigma[[i]][,,k] <- S
        }
        
        # Determine selected features based on mean differences and covariance structure
        mean_diff <- abs(Mu[[i]][k,])
        cov_effect <- colSums(abs(Sigma[[i]][,,k]))
        selectZ_results[[i]][k,] <- (mean_diff > 0.001) | (cov_effect > 0.001)
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

Mstep_XtoY <- function(Y, r, K, N, family, CoY) {
  family <- to_parallel_family(family)


  # if 2 omics layers
  if(length(K) == 2) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(rowSums(lastInd(r,i)), colSums(lastInd(r,i)))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1)]

    if(family == "gaussian") {
      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit))
        Set0 <- cbind(Set0, CoY)
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family = gaussian)
        beta_f <- coef(fit)
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
      }}

    if(family == "binomial") {
      if(is.null(CoY)) {
        # Add try-catch with fallback
        fit <- try({
          glm_fit <- glm(Y ~ r_fit, family = "binomial", control = list(maxit = 25))
          if (!glm_fit$converged) {
            # If not converged, try simpler model
            glm_fit <- glm(Y ~ 1, family = "binomial")
          }
          glm_fit
        }, silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # If error, use intercept-only model
          fit <- glm(Y ~ 1, family = "binomial")
        }
        mu <- as.numeric(coef(fit))
      } else {
        Set0 <- as.data.frame(cbind(Y, r_fit, CoY))
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        
        # Add try-catch with fallback
        fit <- try({
          formula <- as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+")))
          glm_fit <- glm(formula, data = Set0, family = "binomial", control = list(maxit = 25))
          if (!glm_fit$converged) {
            # If not converged, try simpler model
            glm_fit <- glm(Y ~ 1, data = Set0, family = "binomial")
          }
          glm_fit
        }, silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # If error, use intercept-only model
          fit <- glm(Y ~ 1, data = Set0, family = "binomial")
        }
        mu <- as.numeric(coef(fit))
      }

      # Add numerical stability to probability calculation
      mu_array <- vec_to_array(K = K, mu = mu)
      # Subtract max for numerical stability
      mu_array_centered <- mu_array - max(mu_array)
      exp_mu <- exp(mu_array_centered)
      p <- exp_mu / sum(exp_mu)
      
      # Ensure probabilities are in valid range
      p[p < 1e-10] <- 1e-10
      p[p > (1 - 1e-10)] <- 1 - 1e-10
      
      # Renormalize after clamping
      p <- p / sum(p)
      
      sd <- NULL
      mu <- p
    }

    if(any(is.na(mu))) {
      na_indices <- which(is.na(mu))
      # Check which layer has NA values
      layer_with_na <- if(all(na_indices <= K[1])) {
        "Z1"
      } else if (all(na_indices <= sum(K[1:2]))) {
        "Z2"
      } else {
        "later layers"
      }
      stop("No cluster structure is defined for ", layer_with_na)
    }
  }


  # if 3 omics layers
  if(length(K) == 3) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(lastInd(r,i), margin = 1),
        marginSums(lastInd(r,i), margin = 2),
        marginSums(lastInd(r,i), margin = 3))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1)]

    if(family == "gaussian") {
      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
      }else{
        Set0 <- as.data.frame(cbind(Y, r_fit))
        Set0 <- cbind(Set0, CoY)
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        fit <- glm(as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+"))), data = Set0, family = gaussian)

        beta_f <- coef(fit)
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
      }}

    if(family == "binomial") {
      if(is.null(CoY)) {
        # Add try-catch with fallback
        fit <- try({
          glm_fit <- glm(Y ~ r_fit, family = "binomial", control = list(maxit = 25))
          if (!glm_fit$converged) {
            # If not converged, try simpler model
            glm_fit <- glm(Y ~ 1, family = "binomial")
          }
          glm_fit
        }, silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # If error, use intercept-only model
          fit <- glm(Y ~ 1, family = "binomial")
        }
        mu <- as.numeric(coef(fit))
      } else {
        Set0 <- as.data.frame(cbind(Y, r_fit, CoY))
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        
        # Add try-catch with fallback
        fit <- try({
          formula <- as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+")))
          glm_fit <- glm(formula, data = Set0, family = "binomial", control = list(maxit = 25))
          if (!glm_fit$converged) {
            # If not converged, try simpler model
            glm_fit <- glm(Y ~ 1, data = Set0, family = "binomial")
          }
          glm_fit
        }, silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # If error, use intercept-only model
          fit <- glm(Y ~ 1, data = Set0, family = "binomial")
        }
        mu <- as.numeric(coef(fit))
      }

      # Add numerical stability to probability calculation
      mu_array <- vec_to_array(K = K, mu = mu)
      # Subtract max for numerical stability
      mu_array_centered <- mu_array - max(mu_array)
      exp_mu <- exp(mu_array_centered)
      p <- exp_mu / sum(exp_mu)
      
      # Ensure probabilities are in valid range
      p[p < 1e-10] <- 1e-10
      p[p > (1 - 1e-10)] <- 1 - 1e-10
      
      # Renormalize after clamping
      p <- p / sum(p)
      
      sd <- NULL
      mu <- p
    }

    if(any(is.na(mu))) {
      na_indices <- which(is.na(mu))
      # Check which layer has NA values
      layer_with_na <- if(all(na_indices <= K[1])) {
        "Z1"
      } else if (all(na_indices <= sum(K[1:2]))) {
        "Z2"
      } else {
        "later layers"
      }
      stop("No cluster structure is defined for ", layer_with_na)
    }
  }


  # if 4 omics layers
  if(length(K) == 4) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(lastInd(r,i), margin = 1),
        marginSums(lastInd(r,i), margin = 2),
        marginSums(lastInd(r,i), margin = 3),
        marginSums(lastInd(r,i), margin = 4))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1)]

    if(family == "gaussian") {
      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
      }else{
        fit <- lm(Y ~ r_fit + CoY)
        beta_f <- coef(fit)
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
      }}

    if(family == "binomial") {
      if(is.null(CoY)) {
        # Add try-catch with fallback
        fit <- try({
          glm_fit <- glm(Y ~ r_fit, family = "binomial", control = list(maxit = 25))
          if (!glm_fit$converged) {
            # If not converged, try simpler model
            glm_fit <- glm(Y ~ 1, family = "binomial")
          }
          glm_fit
        }, silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # If error, use intercept-only model
          fit <- glm(Y ~ 1, family = "binomial")
        }
        mu <- as.numeric(coef(fit))
      } else {
        Set0 <- as.data.frame(cbind(Y, r_fit, CoY))
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        
        # Add try-catch with fallback
        fit <- try({
          formula <- as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+")))
          glm_fit <- glm(formula, data = Set0, family = "binomial", control = list(maxit = 25))
          if (!glm_fit$converged) {
            # If not converged, try simpler model
            glm_fit <- glm(Y ~ 1, data = Set0, family = "binomial")
          }
          glm_fit
        }, silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # If error, use intercept-only model
          fit <- glm(Y ~ 1, data = Set0, family = "binomial")
        }
        mu <- as.numeric(coef(fit))
      }

      # Add numerical stability to probability calculation
      mu_array <- vec_to_array(K = K, mu = mu)
      # Subtract max for numerical stability
      mu_array_centered <- mu_array - max(mu_array)
      exp_mu <- exp(mu_array_centered)
      p <- exp_mu / sum(exp_mu)
      
      # Ensure probabilities are in valid range
      p[p < 1e-10] <- 1e-10
      p[p > (1 - 1e-10)] <- 1 - 1e-10
      
      # Renormalize after clamping
      p <- p / sum(p)
      
      sd <- NULL
      mu <- p
    }

    if(any(is.na(mu))) {
      na_indices <- which(is.na(mu))
      # Check which layer has NA values
      layer_with_na <- if(all(na_indices <= K[1])) {
        "Z1"
      } else if (all(na_indices <= sum(K[1:2]))) {
        "Z2"
      } else {
        "later layers"
      }
      stop("No cluster structure is defined for ", layer_with_na)
    }
  }


  # if 5 omics layers
  if(length(K) == 5) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(lastInd(r,i), margin = 1),
        marginSums(lastInd(r,i), margin = 2),
        marginSums(lastInd(r,i), margin = 3),
        marginSums(lastInd(r,i), margin = 4),
        marginSums(lastInd(r,i), margin = 5))
    }))
    r_fit <- r_matrix[, -c(1,
                           K[1] + 1,
                           K[1] + K[2] + 1,
                           K[1] + K[2] + K[3] + 1,
                           K[1] + K[2] + K[3] + K[4] + 1)]

    if(family == "gaussian") {
      if(is.null(CoY)) {
        fit <- lm(Y ~ r_fit)
        mu <- as.numeric(coef(fit))
        sd <- sd(resid(fit))
      }else{
        fit <- lm(Y ~ r_fit + CoY)
        beta_f <- summary(fit)$coefficients[, 1]
        mu <- as.numeric(beta_f)
        sd <- sd(resid(fit))
      }}

    if(family == "binomial") {
      if(is.null(CoY)) {
        # Add try-catch with fallback
        fit <- try({
          glm_fit <- glm(Y ~ r_fit, family = "binomial", control = list(maxit = 25))
          if (!glm_fit$converged) {
            # If not converged, try simpler model
            glm_fit <- glm(Y ~ 1, family = "binomial")
          }
          glm_fit
        }, silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # If error, use intercept-only model
          fit <- glm(Y ~ 1, family = "binomial")
        }
        mu <- as.numeric(coef(fit))
      } else {
        Set0 <- as.data.frame(cbind(Y, r_fit, CoY))
        colnames(Set0) <- c("Y", paste0("LC", 1:ncol(r_fit)), colnames(CoY))
        
        # Add try-catch with fallback
        fit <- try({
          formula <- as.formula(paste("Y~", paste(colnames(Set0)[-1], collapse = "+")))
          glm_fit <- glm(formula, data = Set0, family = "binomial", control = list(maxit = 25))
          if (!glm_fit$converged) {
            # If not converged, try simpler model
            glm_fit <- glm(Y ~ 1, data = Set0, family = "binomial")
          }
          glm_fit
        }, silent = TRUE)
        
        if (inherits(fit, "try-error")) {
          # If error, use intercept-only model
          fit <- glm(Y ~ 1, data = Set0, family = "binomial")
        }
        mu <- as.numeric(coef(fit))
      }

      # Add numerical stability to probability calculation
      mu_array <- vec_to_array(K = K, mu = mu)
      # Subtract max for numerical stability
      mu_array_centered <- mu_array - max(mu_array)
      exp_mu <- exp(mu_array_centered)
      p <- exp_mu / sum(exp_mu)
      
      # Ensure probabilities are in valid range
      p[p < 1e-10] <- 1e-10
      p[p > (1 - 1e-10)] <- 1 - 1e-10
      
      # Renormalize after clamping
      p <- p / sum(p)
      
      sd <- NULL
      mu <- p
    }

    if(any(is.na(mu))) {
      na_indices <- which(is.na(mu))
      # Check which layer has NA values
      layer_with_na <- if(all(na_indices <= K[1])) {
        "Z1"
      } else if (all(na_indices <= sum(K[1:2]))) {
        "Z2"
      } else {
        "later layers"
      }
      stop("No cluster structure is defined for ", layer_with_na)
    }
  }



  Delta <- list(mu = mu,
                sd = sd,
                K = K)
  return(list(fit = fit,
              Gamma = Delta))
}

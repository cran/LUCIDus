#' @title Check missing patterns in omics data
#' @param Z A data matrix or list of matrices for parallel model
#' @param lucid_model Model type ("early" or "parallel")
#' @return List containing missing data analysis
#' @export
check_na <- function(Z, lucid_model = c("early", "parallel")) {
  lucid_model <- match.arg(lucid_model)
  
  if (lucid_model == "early") {
    # Original functionality for early model
    missing_analysis <- analyze_missing_pattern(Z)
    
    N <- nrow(Z)
    M <- ncol(Z)
    index <- !is.na(Z)
    obs_na <- rowSums(!index)
    indicator_na <- sapply(1:N, function(i) {
      return(ifelse(obs_na[i] == 0, 1,    # Pattern 1: No missing
                    ifelse(obs_na[i] == M, 3,    # Pattern 3: All missing
                          2)))                    # Pattern 2: Some missing
    })
    impute_flag <- sum(indicator_na == 2) != 0
    
    if (length(missing_analysis$high_miss_cols) > 0) {
      warning("High missingness (>50%) in columns: ", 
              paste(missing_analysis$high_miss_cols, collapse=", "))
    }
    
    return(list(
      index = index,
      indicator_na = indicator_na,
      impute_flag = impute_flag,
      missing_analysis = missing_analysis
    ))
    
  } else {
    # Enhanced functionality for parallel model
    if (!is.list(Z)) {
      stop("For parallel model, Z must be a list of matrices")
    }
    
    # Check consistent number of observations
    N <- nrow(Z[[1]])
    if (!all(sapply(Z, nrow) == N)) {
      stop("All omics layers must have the same number of observations")
    }
    
    # Analyze each layer
    n_layers <- length(Z)
    layer_analyses <- vector("list", n_layers)
    indices <- vector("list", n_layers)
    indicators <- vector("list", n_layers)
    impute_flags <- logical(n_layers)
    
    for (i in 1:n_layers) {
      # Analyze this layer
      layer_analysis <- analyze_missing_pattern(Z[[i]])
      layer_analyses[[i]] <- layer_analysis
      
      # Calculate missing patterns
      index_i <- !is.na(Z[[i]])
      obs_na <- rowSums(!index_i)
      M <- ncol(Z[[i]])
      
      indicator_i <- sapply(1:N, function(j) {
        return(ifelse(obs_na[j] == 0, 1,    # Pattern 1: No missing
                      ifelse(obs_na[j] == M, 3,    # Pattern 3: All missing
                            2)))                    # Pattern 2: Some missing
      })
      
      indices[[i]] <- index_i
      indicators[[i]] <- indicator_i
      impute_flags[i] <- sum(indicator_i == 2) != 0
      
      # Warning for high missingness
      if (length(layer_analysis$high_miss_cols) > 0) {
        warning(sprintf("Layer %d: High missingness (>50%%) in columns: %s", 
                i, paste(layer_analysis$high_miss_cols, collapse=", ")))
      }
    }
    
    # Cross-layer analysis
    cross_layer_patterns <- table(do.call(paste, lapply(indicators, as.character)))
    
    return(list(
      index = indices,
      indicator_na = indicators,
      impute_flag = impute_flags,
      layer_analyses = layer_analyses,
      cross_layer_summary = list(
        n_layers = n_layers,
        n_observations = N,
        features_per_layer = sapply(Z, ncol),
        missing_pattern_counts = cross_layer_patterns,
        total_missing_prop = mean(sapply(layer_analyses, 
                                       function(x) x$total_missing))
      )
    ))
  }
}

#' @title Summarize missing-data patterns from check_na output
#' @param na_pattern Output list from check_na
#' @param lucid_model Model type ("early" or "parallel")
#' @return Missing-data summary list
summarize_missing_stats <- function(na_pattern, lucid_model = c("early", "parallel")) {
  lucid_model <- match.arg(lucid_model)

  summarize_one <- function(index_mat, indicator_vec) {
    n_rows <- nrow(index_mat)
    n_features <- ncol(index_mat)
    listwise_rows <- sum(indicator_vec == 3)
    sporadic_rows <- sum(indicator_vec == 2)
    complete_rows <- sum(indicator_vec == 1)
    total_missing_cells <- sum(!index_mat)
    listwise_missing_cells <- listwise_rows * n_features
    sporadic_missing_cells <- total_missing_cells - listwise_missing_cells

    list(
      n_rows = n_rows,
      n_features = n_features,
      complete_rows = complete_rows,
      listwise_rows = listwise_rows,
      sporadic_rows = sporadic_rows,
      prop_listwise_rows = listwise_rows / n_rows,
      prop_sporadic_rows = sporadic_rows / n_rows,
      total_missing_cells = total_missing_cells,
      sporadic_missing_cells = sporadic_missing_cells,
      prop_total_missing_cells = total_missing_cells / (n_rows * n_features)
    )
  }

  if (lucid_model == "early") {
    return(summarize_one(na_pattern$index, na_pattern$indicator_na))
  }

  # parallel
  n_layers <- length(na_pattern$index)
  layer_summary <- vector("list", n_layers)
  for (i in seq_len(n_layers)) {
    layer_summary[[i]] <- summarize_one(na_pattern$index[[i]], na_pattern$indicator_na[[i]])
  }

  layer_df <- data.frame(
    layer = seq_len(n_layers),
    n_rows = sapply(layer_summary, `[[`, "n_rows"),
    n_features = sapply(layer_summary, `[[`, "n_features"),
    complete_rows = sapply(layer_summary, `[[`, "complete_rows"),
    listwise_rows = sapply(layer_summary, `[[`, "listwise_rows"),
    sporadic_rows = sapply(layer_summary, `[[`, "sporadic_rows"),
    prop_listwise_rows = sapply(layer_summary, `[[`, "prop_listwise_rows"),
    prop_sporadic_rows = sapply(layer_summary, `[[`, "prop_sporadic_rows"),
    total_missing_cells = sapply(layer_summary, `[[`, "total_missing_cells"),
    sporadic_missing_cells = sapply(layer_summary, `[[`, "sporadic_missing_cells"),
    prop_total_missing_cells = sapply(layer_summary, `[[`, "prop_total_missing_cells")
  )

  return(list(
    n_layers = n_layers,
    layer_summary = layer_df
  ))
}

#' @title I-step of LUCID
#' @description Impute missing data in Z by maximizing the likelihood given fixed parameters
#' @param Z Matrix or list of matrices for omics data
#' @param p Probability matrix/vector or list for parallel model
#' @param mu Mean matrix or list of matrices
#' @param sigma Covariance array or list of arrays
#' @param index Missing value indicators
#' @param lucid_model Model type
#' @return Imputed dataset
Istep_Z <- function(Z, p, mu, sigma, index, lucid_model) {
  if (lucid_model == "early") {
    # Early model code
    N <- nrow(Z)
    mu <- t(mu)  # Important: Transpose mu for early model
    
    # Get missing patterns
    na_pattern <- check_na(Z, lucid_model = "early")
    
    # Process each observation
    Z_fill <- t(sapply(1:N, function(i) {
      curr_p <- if (is.matrix(p)) p[i,] else p
      fill_data(obs = Z[i,], mu = mu, sigma = sigma, 
                p = curr_p, index = index[i,], 
                lucid_model = lucid_model)
    }))
    
    return(Z_fill)
    
  } else if (lucid_model == "parallel") {
    # For parallel model, we can receive either a single layer or all layers
    if (is.list(Z) && !is.matrix(Z)) {
      # Multiple layers case - validate all inputs are lists
      if (!is.list(mu) || !is.list(sigma) || !is.list(index)) {
        stop("For parallel model with multiple layers, all inputs must be lists")
      }
      
      n_layers <- length(Z)
      if (!all(sapply(list(mu, sigma, index), length) == n_layers)) {
        stop("Inconsistent number of layers across inputs")
      }
      
      # Validate p structure
      if (is.list(p)) {
        if (length(p) != n_layers) {
          stop("Number of probability sets must match number of layers")
        }
      } else {
        # If p is not a list, replicate it for each layer
        p <- replicate(n_layers, p, simplify = FALSE)
      }
      
      # Process each layer
      Z_fill <- vector("list", n_layers)
      for (i in 1:n_layers) {
        N_i <- nrow(Z[[i]])
        M_i <- ncol(Z[[i]])
        K_i <- nrow(mu[[i]])  # Changed from ncol to nrow for parallel model
        
        # Validate sigma dimensions for this layer
        if (!identical(dim(sigma[[i]]), c(M_i, M_i, K_i))) {
          stop(sprintf("Layer %d: Invalid sigma dimensions. Expected c(%d,%d,%d), got c(%s)", 
                      i, M_i, M_i, K_i, paste(dim(sigma[[i]]), collapse=",")))
        }
        
        # Get missing patterns for this layer
        na_pattern <- check_na(Z[[i]], lucid_model = "early")
        
        # Process observations based on missing pattern
        Z_fill[[i]] <- t(sapply(1:N_i, function(j) {
          pattern <- na_pattern$indicator_na[j]
          curr_p <- if (is.matrix(p[[i]])) p[[i]][j,] else p[[i]]
          
          if (pattern == 1) {
            return(Z[[i]][j,])
          } else if (pattern == 3) {
            return(rep(NA, M_i))
          } else {
            fill_data(obs = Z[[i]][j,], 
                     mu = mu[[i]], 
                     sigma = sigma[[i]], 
                     p = curr_p, 
                     index = na_pattern$index[j,], 
                     lucid_model = lucid_model)
          }
        }))
      }
      
      return(Z_fill)
      
    } else {
      # Single layer case - treat as matrix
      N <- nrow(Z)
      M <- ncol(Z)
      K <- nrow(mu)  # Changed from ncol to nrow for parallel model
      
      # Validate sigma dimensions for parallel model
      if (!identical(dim(sigma), c(M, M, K))) {
        stop(sprintf("Invalid sigma dimensions. Expected c(%d,%d,%d), got c(%s)", 
                    M, M, K, paste(dim(sigma), collapse=",")))
      }
      
      # Get missing patterns
      na_pattern <- check_na(Z, lucid_model = "early")
      
      # Process observations
      Z_fill <- t(sapply(1:N, function(i) {
        pattern <- na_pattern$indicator_na[i]
        curr_p <- if (is.matrix(p)) p[i,] else p
        
        if (pattern == 1) {
          return(Z[i,])
        } else if (pattern == 3) {
          return(rep(NA, M))
        } else {
          fill_data(obs = Z[i,], 
                   mu = mu, 
                   sigma = sigma, 
                   p = curr_p, 
                   index = na_pattern$index[i,], 
                   lucid_model = lucid_model)
        }
      }))
      
      return(Z_fill)
    }
  }
}

#' @title Impute missing data by optimizing the likelihood function
#' @param obs a vector of length M
#' @param mu a matrix of size K x M for parallel model, M x K for early model
#' @param sigma a matrix of size M x M x K
#' @param p a vector of length K
#' @param index a vector of length M, indicating whether a value is missing
#' @param lucid_model Specifying LUCID model
#' @return an observation with updated imputed values
#' @export
fill_data <- function(obs, mu, sigma, p, index, lucid_model) {
  # Input validation
  if (is.null(p) || length(p) == 0) {
    stop("Cluster probabilities cannot be NULL or empty")
  }
  
  # Get dimensions differently based on model type
  if (lucid_model == "parallel") {
    K <- nrow(mu)  # Number of clusters (rows in parallel model)
    M <- ncol(mu)  # Number of features (columns in parallel model)
  } else {
    K <- ncol(mu)  # Number of clusters (columns in early model)
    M <- nrow(mu)  # Number of features (rows in early model)
  }
  
  # Validate dimensions
  if (length(obs) != M) {
    stop(sprintf("Length of obs (%d) must match number of features (%d)", 
                length(obs), M))
  }
  if (length(p) != K) {
    stop(sprintf("Length of probability vector (%d) must match number of clusters (%d)", 
                length(p), K))
  }
  if (length(index) != M) {
    stop(sprintf("Length of index (%d) must match number of features (%d)", 
                length(index), M))
  }
  if (!identical(dim(sigma), c(M, M, K))) {
    stop(sprintf("Invalid sigma dimensions. Expected c(%d,%d,%d), got c(%s)", 
                M, M, K, paste(dim(sigma), collapse=",")))
  }
  
  # Normalize probabilities if needed
  if (abs(sum(p) - 1) > 0.01) {
    warning("Cluster probabilities do not sum to 1, normalizing")
    p <- p / sum(p)
  }
  
  # If no missing values, return original
  if (!any(!index)) {
    return(obs)
  }
  
  # For each missing value, compute weighted average across clusters
  miss_idx <- which(!index)
  for (i in miss_idx) {
    # Extract relevant means and variances for this position
    if (lucid_model == "parallel") {
      cluster_means <- mu[,i]  # For parallel model, means are in columns
    } else {
      cluster_means <- mu[i,]  # For early model, means are in rows
    }
    cluster_vars <- sapply(1:K, function(k) sigma[i,i,k])
    
    # Get observed values for this feature (if any)
    obs_idx <- which(index)
    if (length(obs_idx) > 0) {
      # Use conditional distribution if we have observed values
      weights <- try({
        # Compute log-likelihoods for each cluster
        log_likes <- sapply(1:K, function(k) {
          if (!is.finite(cluster_vars[k]) || cluster_vars[k] <= 0) {
            return(-Inf)
          }
          # Use full covariance for likelihood calculation
          if (lucid_model == "parallel") {
            obs_mean <- mu[k, obs_idx]  # For parallel model
          } else {
            obs_mean <- mu[obs_idx, k]  # For early model
          }
          obs_sigma <- sigma[obs_idx, obs_idx, k]
          obs_val <- obs[obs_idx]
          
          tryCatch({
            dmvnorm(obs_val, mean = obs_mean, sigma = obs_sigma, log = TRUE)
          }, error = function(e) {
            warning("Error in mvnorm calculation, using diagonal approximation")
            sum(dnorm(obs_val, mean = obs_mean, 
                     sd = sqrt(diag(obs_sigma)), log = TRUE))
          })
        })
        
        # Combine with prior probabilities
        log_post <- log_likes + log(p)
        
        # Normalize to get posterior probabilities
        exp(log_post - max(log_post)) / sum(exp(log_post - max(log_post)))
      }, silent = TRUE)
    } else {
      # No observed values, use prior probabilities
      weights <- p
    }
    
    # If computation failed, use prior probabilities
    if (inherits(weights, "try-error")) {
      weights <- p
    }
    
    # Impute using weighted average
    obs[i] <- sum(cluster_means * weights)
  }
  
  return(obs)
}

# Calculate the first half of the imputed values with enhanced stability
fill_data_help1 <- function(obs, B, mu, alpha, sigma_inv, P) {
  K <- ncol(mu)
  l <- length(B)
  res <- matrix(rep(0, l * l), nrow = l)
  
  # Accumulate sum with stability check
  for(j in 1:K) {
    contrib <- alpha[j] * P[j] * sigma_inv[,,j][B, B]
    if (any(is.infinite(contrib))) {
      warning("Infinite values detected in matrix calculation")
      contrib[is.infinite(contrib)] <- sign(contrib[is.infinite(contrib)]) * 1e10
    }
    res <- res + contrib
  }
  
  # Safe matrix inversion
  res_inv <- safe_solve(res)
  if (is.null(res_inv)) {
    warning("Matrix inversion failed in fill_data_help1, using regularized version")
    res <- res + diag(1e-6, nrow(res))
    res_inv <- solve(res)
  }
  
  return(res_inv)
}

# Calculate the second half of the imputed values with enhanced stability
fill_data_help2 <- function(obs, A, B, mu, alpha, sigma_inv, P) {
  K <- ncol(mu)
  l <- length(B)
  res <- rep(0, l)
  
  for(j in 1:K) {
    s1 <- sigma_inv[,,j][B, A, drop = FALSE]
    s2 <- sigma_inv[,,j][B, B, drop = FALSE]
    mu1 <- mu[A, j, drop = FALSE]
    mu2 <- mu[B, j, drop = FALSE]
    
    # Safe matrix multiplication with checks
    term1 <- tryCatch({
      s1 %*% mu1
    }, error = function(e) {
      warning("Error in matrix multiplication 1: ", e$message)
      matrix(0, nrow = nrow(s1), ncol = 1)
    })
    
    term2 <- tryCatch({
      s2 %*% mu2
    }, error = function(e) {
      warning("Error in matrix multiplication 2: ", e$message)
      matrix(0, nrow = nrow(s2), ncol = 1)
    })
    
    term3 <- tryCatch({
      s1 %*% obs[A]
    }, error = function(e) {
      warning("Error in matrix multiplication 3: ", e$message)
      matrix(0, nrow = nrow(s1), ncol = 1)
    })
    
    xx <- term1 + term2 - term3
    
    # Check for numerical stability
    if (any(is.infinite(xx))) {
      warning("Infinite values in calculation, using bounded values")
      xx[is.infinite(xx)] <- sign(xx[is.infinite(xx)]) * 1e10
    }
    
    res <- res + alpha[j] * P[j] * xx
  }
  
  return(res)
}

# impute missing values in Z by LOD
fill_data_lod <- function(Z_vec) {
  na_ind <- is.na(Z_vec)
  if(any(na_ind)) {
    if(!all(na_ind)) {
      lod <- min(Z_vec, na.rm = TRUE)
      Z_vec[na_ind] <- lod / sqrt(2)
    }
  }
  return(Z_vec)
}

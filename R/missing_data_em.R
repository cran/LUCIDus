#' @title Classify each subject's omics missingness pattern
#' @description
#' Assigns every subject to one of the three missingness categories the LUCID EM
#' algorithm branches on, and reports whether imputation is required at all.
#' The categories follow the incomplete-omics extension of LUCID:
#' \describe{
#'   \item{complete (code 1)}{Every feature observed. Contributes the ordinary
#'     complete-data term to the likelihood.}
#'   \item{sporadic (code 2)}{Some but not all features observed. These are the
#'     subjects that require imputation: the missing coordinates are integrated
#'     out against the fitted cluster model in the I-step.}
#'   \item{listwise (code 3)}{No feature observed for this layer. The omics term
#'     drops out of that subject's likelihood entirely, so the subject still
#'     informs the exposure and outcome models but needs no imputation.}
#' }
#'
#' The distinction matters for cost as well as correctness: \code{impute_flag} is
#' \code{TRUE} only when at least one \emph{sporadic} subject exists, and a
#' dataset whose missingness is purely listwise is fitted without any imputation
#' step.
#'
#' A warning is issued for any feature more than half missing, per layer.
#'
#' @param Z For \code{lucid_model = "early"}, an N by M omics matrix. For
#'   \code{"parallel"}, a list of such matrices, one per layer, all with the same
#'   number of rows; anything else is an error.
#' @param lucid_model Either \code{"early"} or \code{"parallel"}. A serial fit
#'   calls this once per stage rather than passing \code{"serial"} here.
#' @return For \code{"early"}, a list with \code{index} (an N by M logical matrix
#'   that is \code{TRUE} where observed), \code{indicator_na} (the length-N
#'   vector of codes 1, 2, 3 above), \code{impute_flag} (a single logical), and
#'   \code{missing_analysis} (the \code{\link{analyze_missing_pattern}} result).
#'
#'   For \code{"parallel"}, \code{index}, \code{indicator_na} and
#'   \code{layer_analyses} are lists with one element per layer,
#'   \code{impute_flag} is a logical vector over layers, and
#'   \code{cross_layer_summary} adds \code{n_layers}, \code{n_observations},
#'   \code{features_per_layer}, \code{missing_pattern_counts} (a table of the
#'   joint across-layer pattern, so that subjects missing an entire layer can be
#'   distinguished from those missing scattered features in several) and
#'   \code{total_missing_prop}, the proportion of missing \emph{cells} pooled
#'   over layers rather than the mean of per-layer rates, which would weight a
#'   one-feature layer as heavily as a fifty-feature one.
#' @seealso \code{\link{analyze_missing_pattern}} for the per-layer detail.
#' @examples
#' Z <- matrix(rnorm(200), nrow = 20)
#' Z[1:2, 1] <- NA   # sporadic
#' Z[20, ] <- NA     # listwise
#' table(check_na(Z, lucid_model = "early")$indicator_na)
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
        # D9: proportion of missing CELLS across all layers.  Averaging the
        # per-layer rates weights a 1-feature layer the same as a 50-feature
        # layer and is not a proportion of anything.
        total_missing_prop = sum(sapply(indices, function(x) sum(!x))) /
                             sum(sapply(indices, length))
      )
    ))
  }
}

#' @title Summarize missing-data patterns from check_na output
#' @description Reduces the per-subject classification from \code{check_na} to
#'   the row and cell counts stored as a fitted model's \code{missing_summary}.
#' @param na_pattern Output list from \code{check_na}.
#' @param lucid_model Model type ("early" or "parallel").
#' @return Missing-data summary list: for "early" a single list of counts and
#'   proportions, for "parallel" one such list per layer.
#' @noRd
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
#' @noRd
Istep_Z <- function(Z, p, mu, sigma, index, lucid_model) {
  if (lucid_model == "early") {
    # Early model code
    N <- nrow(Z)
    mu <- t(mu)  # Important: Transpose mu for early model
    
    # Process each observation
    Z_fill <- t(sapply(1:N, function(i) {
      curr_p <- if (is.matrix(p)) p[i,] else p
      if (all(index[i, ])) return(Z[i, ])
      if (!any(index[i, ])) return(Z[i, ])
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
        
        # Process observations based on missing pattern
        Z_fill[[i]] <- t(sapply(1:N_i, function(j) {
          curr_p <- if (is.matrix(p[[i]])) p[[i]][j,] else p[[i]]
          
          if (all(index[[i]][j, ])) {
            return(Z[[i]][j,])
          } else if (!any(index[[i]][j, ])) {
            return(Z[[i]][j,])
          } else {
            fill_data(obs = Z[[i]][j,], 
                     mu = mu[[i]], 
                     sigma = sigma[[i]], 
                     p = curr_p, 
                     index = index[[i]][j,], 
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
      
      # Process observations
      Z_fill <- t(sapply(1:N, function(i) {
        curr_p <- if (is.matrix(p)) p[i,] else p
        
        if (all(index[i, ])) {
          return(Z[i,])
        } else if (!any(index[i, ])) {
          return(Z[i,])
        } else {
          fill_data(obs = Z[i,], 
                   mu = mu, 
                   sigma = sigma, 
                   p = curr_p, 
                   index = index[i,], 
                   lucid_model = lucid_model)
        }
      }))
      
      return(Z_fill)
    }
  }
}

#' @title Impute one subject's missing omics values under the cluster model
#' @description
#' The I-step for a single subject. Given the current cluster-specific means and
#' covariances and the subject's cluster probabilities, the missing coordinates
#' are set to the values that maximize the mixture log-likelihood, holding the
#' observed coordinates and all parameters fixed.
#'
#' Writing \eqn{A} for the observed coordinates, \eqn{B} for the missing ones and
#' \eqn{\Omega^{k}} for cluster \eqn{k}'s precision matrix, the stationary point
#' of the weighted log-density in \eqn{x_B} is the solution of
#' \deqn{\Big(\sum_k w_k \Omega^{k}_{BB}\Big) x_B = \sum_k w_k \Big(
#'   \Omega^{k}_{BB} \mu^{k}_{B} + \Omega^{k}_{BA} \mu^{k}_{A}
#'   - \Omega^{k}_{BA} x_A \Big),}
#' which is what this function solves. It is a precision-weighted combination of
#' the clusters, not an average of their individual conditional means: a cluster
#' with a sharply concentrated covariance pulls the imputed value harder than a
#' diffuse one. For \eqn{K = 1} it reduces exactly to the Gaussian conditional
#' mean \eqn{\mu_B - \Omega_{BB}^{-1}\Omega_{BA}(x_A - \mu_A)}.
#'
#' The weights \eqn{w_k} are not \code{p} itself. They are \code{p} reweighted by
#' each cluster's density at the subject's current values and renormalised, so
#' that clusters the subject's observed coordinates are incompatible with are
#' discounted; if no cluster yields a finite density, \code{p} is used unchanged.
#' Any missing entry still \code{NA} on entry is first given the \code{p}-weighted
#' mean of the cluster means, so the density evaluation has a complete vector to
#' work on.
#'
#' Conditioning on the observed coordinates -- rather than substituting a cluster
#' mean outright -- is what keeps the imputation consistent with the correlation
#' structure the model has estimated, and is the step that makes the
#' missing-at-random assumption of the incomplete-omics extension operative.
#'
#' A subject with nothing observed, or nothing missing, is returned unchanged:
#' there is no conditional information in the first case and nothing to do in the
#' second. Fully unobserved subjects are handled listwise instead, see
#' \code{\link{check_na}}.
#'
#' @param obs A length-M vector holding one subject's omics values, missing
#'   entries included.
#' @param mu Cluster means: a K by M matrix for the parallel model, M by K for
#'   the early model.
#' @param sigma An M by M by K array of cluster covariances.
#' @param p A length-K vector of that subject's posterior cluster
#'   probabilities, used as the averaging weights.
#' @param index A length-M logical vector, \code{TRUE} where the corresponding
#'   entry of \code{obs} is observed.
#' @param lucid_model Specifying LUCID model, which fixes the orientation of
#'   \code{mu}.
#' @return The length-M vector \code{obs} with its missing entries replaced by
#'   their conditional expectations; observed entries are unchanged.
#' @seealso \code{\link{safe_impute}} for the cruder single-value filling used
#'   to initialize the algorithm before any cluster model exists.
#' @noRd
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
  
  if (any(!is.finite(p)) || sum(p) <= 0) stop("Cluster probabilities must be finite and positive")
  p <- p / sum(p)
  
  A <- which(index)
  B <- which(!index)
  if (!length(B) || !length(A)) return(obs)

  mu_km <- if (lucid_model == "parallel") as.matrix(mu) else t(as.matrix(mu))
  current <- as.numeric(obs)
  if (anyNA(current[B])) current[B] <- as.numeric(p %*% mu_km[, B, drop = FALSE])

  precision <- lapply(seq_len(K), function(k) {
    safe_solve(check_and_stabilize_sigma(sigma[, , k]))
  })
  log_density <- vapply(seq_len(K), function(k) {
    value <- try(mclust::dmvnorm(matrix(current, nrow = 1L), mean = mu_km[k, ],
                                 sigma = sigma[, , k], log = TRUE), silent = TRUE)
    if (inherits(value, "try-error") || !is.finite(value)) -Inf else value
  }, numeric(1))
  log_weight <- log(p) + log_density
  weights <- if (all(!is.finite(log_weight))) p else {
    exp(log_weight - safe_log_sum_exp(log_weight))
  }

  lhs <- matrix(0, length(B), length(B))
  rhs <- numeric(length(B))
  for (k in seq_len(K)) {
    omega_bb <- precision[[k]][B, B, drop = FALSE]
    omega_ba <- precision[[k]][B, A, drop = FALSE]
    lhs <- lhs + weights[k] * omega_bb
    rhs <- rhs + weights[k] * as.numeric(
      omega_ba %*% mu_km[k, A] + omega_bb %*% mu_km[k, B] -
        omega_ba %*% current[A]
    )
  }
  current[B] <- as.numeric(safe_solve(lhs) %*% rhs)
  current
}

#' Impute missing omics values by limit-of-detection substitution
#'
#' Replaces \code{NA} entries in one feature with \eqn{\text{LOD}/\sqrt{2}},
#' where LOD is taken as the observed minimum for that feature. A
#' fully-missing feature is left unchanged (\code{NA} restored, not a
#' fabricated LOD).
#'
#' @param Z_vec One omics feature's values, possibly containing \code{NA}.
#' @return \code{Z_vec} with \code{NA} entries replaced, unless every entry
#'   is missing.
#' @noRd
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

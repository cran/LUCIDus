# Enhanced utilities for missing data handling

#' @title Analyze missing data patterns in detail
#' @param Z Matrix of data to analyze for missing patterns
#' @return List containing detailed missing data analysis
#' @export
analyze_missing_pattern <- function(Z) {
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  
  # Calculate missingness by column
  col_miss <- colMeans(is.na(Z))
  high_miss_cols <- which(col_miss > 0.5)
  
  # Calculate missingness by row
  row_miss <- rowMeans(is.na(Z))
  high_miss_rows <- which(row_miss > 0.5)
  
  # Pattern analysis
  n_obs <- nrow(Z)
  n_vars <- ncol(Z)
  complete_cases <- sum(complete.cases(Z))
  
  # Missing patterns
  patterns <- unique(is.na(Z))
  n_patterns <- nrow(patterns)
  
  return(list(
    col_missingness = col_miss,
    row_missingness = row_miss,
    high_miss_cols = high_miss_cols,
    high_miss_rows = high_miss_rows,
    n_complete = complete_cases,
    n_patterns = n_patterns,
    total_missing = sum(is.na(Z)) / (n_obs * n_vars)
  ))
}

#' @title Check quality of imputed data
#' @param original Original data matrix with missing values
#' @param imputed Imputed data matrix
#' @return List containing imputation quality metrics and overall assessment
#' @export
check_imputation_quality <- function(original, imputed) {
  if (!is.matrix(original)) original <- as.matrix(original)
  if (!is.matrix(imputed)) imputed <- as.matrix(imputed)
  
  # Check dimensions match
  if (any(dim(original) != dim(imputed))) {
    warning("Dimensions of original and imputed data do not match")
    return(list(
      is_valid = FALSE,
      mean_diff = NA,
      sd_ratio = NA,
      warning = "Dimension mismatch"
    ))
  }
  
  # Initialize aggregated metrics
  all_mean_diffs <- numeric(0)
  all_sd_ratios <- numeric(0)
  
  # Analyze each column
  for (j in 1:ncol(original)) {
    # Get non-missing values in original data
    orig_complete <- original[!is.na(original[,j]), j]
    imp_values <- imputed[is.na(original[,j]), j]  # Only look at imputed values
    
    if (length(orig_complete) > 0 && length(imp_values) > 0) {
      # Calculate statistics safely
      orig_mean <- mean(orig_complete, na.rm = TRUE)
      orig_sd <- sd(orig_complete, na.rm = TRUE)
      imp_mean <- mean(imp_values, na.rm = TRUE)
      imp_sd <- sd(imp_values, na.rm = TRUE)
      
      # Handle zero/NA standard deviations
      if (is.na(orig_sd) || orig_sd == 0) {
        if (is.na(imp_sd) || imp_sd == 0) {
          # Both have no variation - this is okay
          mean_diff <- if (is.na(orig_mean) || is.na(imp_mean)) NA else (imp_mean - orig_mean)
          sd_ratio <- 1
        } else {
          # Original has no variation but imputed does - potentially problematic
          mean_diff <- NA
          sd_ratio <- Inf
        }
      } else {
        # Normal case - both have variation
        mean_diff <- (imp_mean - orig_mean) / orig_sd
        sd_ratio <- if (is.na(imp_sd)) NA else (imp_sd / orig_sd)
      }
      
      # Store non-NA values
      if (!is.na(mean_diff)) all_mean_diffs <- c(all_mean_diffs, mean_diff)
      if (!is.na(sd_ratio)) all_sd_ratios <- c(all_sd_ratios, sd_ratio)
    }
  }
  
  # Compute overall metrics
  mean_diff <- if (length(all_mean_diffs) > 0) mean(abs(all_mean_diffs)) else NA
  sd_ratio <- if (length(all_sd_ratios) > 0) mean(all_sd_ratios) else NA
  
  # Final quality assessment
  is_valid <- TRUE
  warning_msg <- NULL
  
  if (is.na(mean_diff) && is.na(sd_ratio)) {
    is_valid <- FALSE
    warning_msg <- "No valid comparisons possible"
  } else {
    if (!is.na(mean_diff) && mean_diff > 2) {
      is_valid <- FALSE
      warning_msg <- "Large differences in means detected"
    }
    if (!is.na(sd_ratio) && (sd_ratio > 3 || sd_ratio < 0.3)) {
      is_valid <- FALSE
      warning_msg <- paste0(warning_msg, if (!is.null(warning_msg)) "; ", "Variance ratios out of acceptable range")
    }
  }
  
  return(list(
    is_valid = is_valid,
    mean_diff = mean_diff,
    sd_ratio = sd_ratio,
    warning = warning_msg
  ))
}

#' Safe imputation for edge cases
#' @param Z Matrix with missing values
#' @param method Imputation method ("mean", "median", "lod")
#' @return Imputed matrix
#' @export
safe_impute <- function(Z, method = c("mean", "median", "lod")) {
  method <- match.arg(method)
  
  # Handle edge cases
  if (all(is.na(Z))) {
    warning("All values are missing, cannot impute")
    return(Z)
  }
  
  # Impute column by column
  Z_imp <- apply(Z, 2, function(x) {
    if (all(is.na(x))) {
      warning("Entire column missing, using global mean/median")
      return(rep(mean(Z, na.rm=TRUE), length(x)))
    }
    
    switch(method,
           "mean" = {
             x[is.na(x)] <- mean(x, na.rm=TRUE)
           },
           "median" = {
             x[is.na(x)] <- median(x, na.rm=TRUE)
           },
           "lod" = {
             lod <- min(x, na.rm=TRUE)
             x[is.na(x)] <- lod/sqrt(2)
           })
    return(x)
  })
  
  return(Z_imp)
} 

# Enhanced utilities for missing data handling

#' @title Describe the missing-data pattern of an omics matrix
#' @description
#' Summarises where missingness sits in an omics layer before a LUCID model is
#' fitted, so that the choice between listwise and sporadic handling can be made
#' from the data rather than assumed. Reports missingness by feature and by
#' subject, flags the features and subjects that are more than half missing, and
#' counts the distinct missingness patterns present.
#'
#' The number of distinct patterns is the diagnostic that matters most for cost:
#' the observed-data likelihood is evaluated once per pattern, so a matrix with
#' few patterns (largely listwise missingness) is far cheaper to fit than one of
#' the same sparsity spread over many patterns.
#'
#' @param Z An N by M omics matrix, or an object coercible to one by
#'   \code{as.matrix}. Missing values are \code{NA}.
#' @return A list with components:
#'   \describe{
#'     \item{col_missingness}{Proportion missing for each of the M features.}
#'     \item{row_missingness}{Proportion missing for each of the N subjects.}
#'     \item{high_miss_cols}{Integer indices of features more than half missing.}
#'     \item{high_miss_rows}{Integer indices of subjects more than half missing.}
#'     \item{n_complete}{Number of subjects with no missing feature.}
#'     \item{n_patterns}{Number of distinct missingness patterns, counting the
#'       complete pattern if any subject is complete.}
#'     \item{total_missing}{Proportion of missing cells over the whole matrix.}
#'   }
#' @seealso \code{\link{check_na}}, which classifies subjects into the complete /
#'   sporadic / listwise categories the EM algorithm branches on.
#' @examples
#' Z <- matrix(rnorm(200), nrow = 20)
#' Z[1:3, 1] <- NA
#' Z[5, ] <- NA
#' analyze_missing_pattern(Z)[c("n_complete", "n_patterns", "total_missing")]
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

#' @title Check whether imputed values are distributionally plausible
#' @description
#' Compares the values that were filled in against the values that were actually
#' observed, feature by feature, and flags an imputation that has shifted the
#' centre or the spread of the data far enough to distort a subsequent LUCID
#' fit. This is a sanity check on the imputation, not a measure of its accuracy:
#' the true values are unknown, so agreement in distribution is all that can be
#' assessed.
#'
#' For each feature, the imputed values are compared with the observed ones
#' through a standardised mean difference,
#' \eqn{(\bar{x}_{imp} - \bar{x}_{obs}) / s_{obs}}, and a spread ratio
#' \eqn{s_{imp} / s_{obs}}. Features with no observed values or no imputed
#' values are skipped. The reported \code{mean_diff} is the mean \emph{absolute}
#' standardised difference across features, and \code{sd_ratio} the mean spread
#' ratio. A feature that is constant where observed contributes a
#' \code{sd_ratio} of 1 if its imputed values are also constant, and \code{Inf}
#' if they are not.
#'
#' The imputation is declared invalid when \code{mean_diff} exceeds 2 (the
#' filled values sit more than two observed standard deviations from the
#' observed centre), or when \code{sd_ratio} falls outside \code{[0.3, 3]} --
#' below that range indicates the near-constant imputation that mean-filling
#' produces at high missingness, which biases cluster covariances towards
#' singularity.
#'
#' @param original The data matrix before imputation, containing \code{NA}.
#' @param imputed The same matrix after imputation, with identical dimensions
#'   and column order. A dimension mismatch is a warning, not an error, and
#'   returns \code{is_valid = FALSE}.
#' @return A list with components:
#'   \describe{
#'     \item{is_valid}{\code{TRUE} if neither threshold was breached and at
#'       least one feature could be compared.}
#'     \item{mean_diff}{Mean absolute standardised mean difference across
#'       comparable features, or \code{NA} if none.}
#'     \item{sd_ratio}{Mean ratio of imputed to observed standard deviation, or
#'       \code{NA} if none.}
#'     \item{warning}{\code{NULL} when valid; otherwise a string naming the
#'       thresholds that were breached.}
#'   }
#' @seealso \code{\link{safe_impute}} for the imputations this is meant to
#'   check.
#' @examples
#' Z <- matrix(rnorm(200), nrow = 20)
#' Z_na <- Z; Z_na[1:5, 1] <- NA
#' check_imputation_quality(Z_na, safe_impute(Z_na, method = "mean"))
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

#' @title Single-value imputation that tolerates fully missing columns
#' @description
#' Fills missing entries feature by feature with a single summary of that
#' feature's observed values. This is a starting point for the EM algorithm, not
#' a substitute for it: LUCID's own E-step imputes missing omics values under
#' the fitted cluster model (an internal EM detail, not part of the public
#' API), and single-value
#' filling here only has to be finite and roughly located.
#'
#' Unlike a bare \code{mean(x, na.rm = TRUE)}, a feature with \emph{no} observed
#' value does not yield \code{NaN}: it falls back to the mean over the whole
#' matrix and warns. Note that this fallback uses the mean whichever
#' \code{method} was requested, since a median or limit of detection is not
#' defined for a feature with nothing observed.
#'
#' @param Z A numeric matrix with missing values coded \code{NA}.
#' @param method The summary used to fill a feature:
#'   \describe{
#'     \item{"mean"}{The feature's observed mean.}
#'     \item{"median"}{The feature's observed median; preferable for a skewed
#'       feature, where the mean is pulled towards the tail.}
#'     \item{"lod"}{The feature's observed minimum divided by \eqn{\sqrt{2}},
#'       the standard substitution for values below an assay's limit of
#'       detection. This is appropriate only for data on the original
#'       measurement scale, where the minimum stands in for the detection limit;
#'       on centred or standardised data the observed minimum is negative and
#'       dividing it by \eqn{\sqrt{2}} moves the filled value \emph{up}, which
#'       is not what the convention intends.}
#'   }
#' @return A matrix of the same dimensions as \code{Z} with no missing values,
#'   unless every value of \code{Z} is missing, in which case \code{Z} is
#'   returned unchanged with a warning.
#' @seealso \code{\link{check_imputation_quality}} to check the result, and the
#'   \code{init_impute} argument of \code{\link{estimate_lucid}} for the
#'   imputation LUCID applies internally.
#' @examples
#' Z <- matrix(rnorm(100), nrow = 10)
#' Z[2:4, 2] <- NA
#' colMeans(is.na(safe_impute(Z, method = "median")))
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

# Utility functions for numerical stability in EM algorithm

#' Check matrix condition and stabilize if needed
#' @param sigma Covariance matrix to check
#' @param threshold Condition number threshold (default 1e10)
#' @param epsilon Small constant for regularization (default 1e-6)
#' @return Stabilized covariance matrix
#' @export
check_and_stabilize_sigma <- function(sigma, threshold = 1e10, epsilon = 1e-6) {
  .lucid_check_and_stabilize_sigma(sigma, threshold = threshold, epsilon = epsilon)
}

#' Safe matrix inversion with stability checks
#' @param matrix Matrix to invert
#' @param tol Tolerance for numerical stabilization
#' @return Inverted matrix or NULL if inversion fails
#' @export
safe_solve <- function(matrix, tol = 1e-6) {
  .lucid_safe_solve(matrix, threshold = 1e10, epsilon = tol)
}

#' Safe log-sum-exp computation
#' @param x Numeric vector
#' @return Numerically stable result
#' @export
safe_log_sum_exp <- function(x) {
  .lucid_safe_log_sum_exp(x)
}

#' Check convergence with both absolute and relative criteria
#' @param old_val Previous value
#' @param new_val New value
#' @param abs_tol Absolute tolerance
#' @param rel_tol Relative tolerance
#' @return Boolean indicating convergence
#' @export
check_convergence <- function(old_val, new_val, abs_tol = 1e-6, rel_tol = 1e-6) {
  # Check for invalid values
  if (!is.finite(old_val) || !is.finite(new_val)) {
    return(FALSE)
  }
  
  # Absolute difference
  abs_diff <- abs(new_val - old_val)
  
  # Relative difference
  rel_diff <- abs_diff / (abs(old_val) + .Machine$double.eps)
  
  # Check both criteria
  abs_conv <- abs_diff < abs_tol
  rel_conv <- rel_diff < rel_tol
  
  return(abs_conv || rel_conv)
}

#' Safe probability normalization
#' @param x Vector of probabilities
#' @return Normalized probabilities
#' @export
safe_normalize <- function(x) {
  .lucid_safe_normalize(x)
} 

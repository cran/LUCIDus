# Utility functions for numerical stability

#' Safe log sum exp operation
#' 
#' @param x Numeric vector or matrix
#' @return Log of sum of exponentials, computed in a numerically stable way
#' @noRd
safe_log_sum_exp <- function(x) {
  .lucid_safe_log_sum_exp(x)
}

#' Safe log operation with minimum value
#' 
#' @param x Numeric value
#' @param min_val Minimum value before taking log (default: 1e-300)
#' @return Safe log value
#' @noRd
safe_log <- function(x, min_val = 1e-300) {
  log(pmax(x, min_val))
}

#' Make matrix positive definite
#' 
#' @param matrix Input matrix
#' @param epsilon Minimum eigenvalue (default: 1e-6)
#' @return Positive definite matrix
#' @noRd
make_positive_definite <- function(matrix, epsilon = 1e-6) {
  if(!is.matrix(matrix)) {
    stop("Input must be a matrix")
  }
  
  # Ensure matrix is symmetric
  matrix <- (matrix + t(matrix)) / 2
  
  # Eigendecomposition
  eigen_decomp <- eigen(matrix, symmetric = TRUE)
  values <- eigen_decomp$values
  vectors <- eigen_decomp$vectors
  
  # Fix negative or small eigenvalues
  values <- pmax(values, epsilon)
  
  # Reconstruct matrix
  result <- vectors %*% diag(values) %*% t(vectors)
  
  # Ensure symmetry (due to numerical precision)
  result <- (result + t(result)) / 2
  
  return(result)
}

#' Check if matrix is positive definite
#' 
#' @param matrix Input matrix
#' @return Logical indicating if matrix is positive definite
#' @noRd
is_positive_definite <- function(matrix) {
  if(!is.matrix(matrix)) return(FALSE)
  if(!isSymmetric(matrix)) return(FALSE)
  
  # Try Cholesky decomposition
  tryCatch({
    chol(matrix)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

#' Safe matrix inversion
#' 
#' @param matrix Input matrix
#' @param tol Tolerance for eigenvalues (default: 1e-6)
#' @return Inverted matrix
#' @noRd
safe_solve <- function(matrix, tol = 1e-6) {
  .lucid_safe_solve(matrix, threshold = 1e10, epsilon = tol)
} 

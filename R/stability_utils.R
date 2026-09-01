# Numerical-stability utilities for the LUCID EM algorithm.
#
# These are the single canonical definitions of the exported `safe_*` /
# `check_*` helpers. Each is a thin, documented wrapper around the
# corresponding internal `.lucid_*` implementation below.

#' Condition a covariance matrix for safe downstream use
#'
#' Repairs a covariance matrix so that the multivariate-normal density and
#' precision computations of the EM algorithm cannot fail on it. Three
#' problems are handled in order: non-finite entries are replaced by the
#' (usually finite) diagonal, with any remaining non-positive variance set to
#' \code{epsilon}; the matrix is symmetrised as \eqn{(\Sigma + \Sigma^{T})/2}
#' to remove asymmetry accumulated by floating-point M-step updates; and a
#' ridge of \code{epsilon} is added to the diagonal if the condition number
#' exceeds \code{threshold} or the smallest eigenvalue falls below
#' \code{epsilon}.
#'
#' Ill-conditioned cluster covariances are routine in LUCID rather than
#' exceptional: high-dimensional omics layers with a small effective sample
#' size in a cluster can drive a covariance towards singularity during EM.
#' Conditioning the matrix keeps the affected cluster in the model instead of
#' aborting the fit.
#'
#' @param sigma A square covariance matrix. A non-matrix or non-square
#'   argument is returned unchanged.
#' @param threshold Condition-number threshold above which the ridge is
#'   applied (default \code{1e10}).
#' @param epsilon Ridge magnitude, and the floor imposed on variances and
#'   eigenvalues (default \code{1e-6}).
#' @return A symmetric, finite, positive-definite matrix of the same dimension
#'   as \code{sigma}, or \code{sigma} unchanged if it is not a square matrix.
#' @seealso \code{\link{safe_solve}}, which applies this conditioning as its
#'   own fallback before inverting.
#' @noRd
check_and_stabilize_sigma <- function(sigma, threshold = 1e10, epsilon = 1e-6) {
  .lucid_check_and_stabilize_sigma(sigma, threshold = threshold, epsilon = epsilon)
}

#' Invert a matrix, falling back through conditioning and pseudo-inversion
#'
#' Attempts \code{\link[base]{solve}} first. If that fails, the matrix is
#' conditioned with \code{\link{check_and_stabilize_sigma}} and inversion is
#' retried. If it fails again, a Moore-Penrose pseudo-inverse is formed from
#' the singular value decomposition, inverting only singular values greater
#' than \code{tol} and setting the reciprocal of the rest to zero.
#'
#' Unlike \code{solve()}, this never signals an error: a matrix that cannot be
#' inverted by any of the three routes yields \code{NULL}, which callers in the
#' EM algorithm test for and treat as a failed cluster update.
#'
#' @param matrix A square matrix to invert. A non-matrix or non-square
#'   argument yields \code{NULL}.
#' @param tol Singular values at or below this value are treated as zero when
#'   the pseudo-inverse fallback is reached; also the ridge magnitude passed to
#'   \code{\link{check_and_stabilize_sigma}} (default \code{1e-6}).
#' @param threshold Condition-number threshold forwarded to
#'   \code{\link{check_and_stabilize_sigma}} (default \code{1e10}).
#' @return The inverse (or pseudo-inverse) of \code{matrix}, or \code{NULL} if
#'   the argument is not a square matrix or all three routes fail.
#' @noRd
safe_solve <- function(matrix, tol = 1e-6, threshold = 1e10) {
  .lucid_safe_solve(matrix, threshold = threshold, epsilon = tol)
}

#' Log-sum-exp, computed without overflow
#'
#' Evaluates \eqn{\log \sum_k \exp(x_k)} by factoring out the maximum, the
#' standard stable form used throughout the E-step to normalise cluster
#' log-joint densities that would otherwise underflow to zero when
#' exponentiated directly.
#'
#' @param x A numeric vector, or a numeric matrix. For a matrix the reduction
#'   is applied across each row, so a vector of row-wise results is returned.
#' @return For a vector, a single numeric value; for a matrix, a numeric vector
#'   of length \code{nrow(x)}. An all-\code{-Inf} input yields \code{-Inf}
#'   rather than \code{NaN}.
#' @noRd
safe_log_sum_exp <- function(x) {
  .lucid_safe_log_sum_exp(x)
}

#' Test EM convergence on absolute or relative change
#'
#' Reports convergence when \strong{either} criterion is met: the absolute
#' change \eqn{|new - old|} is below \code{abs_tol}, or the relative change
#' \eqn{|new - old| / |old|} is below \code{rel_tol}. Either alone is
#' unreliable across the scales a LUCID log-likelihood can take -- an absolute
#' tolerance is unreachable for a log-likelihood of large magnitude, and a
#' relative tolerance is unreachable near zero -- so satisfying one is taken as
#' sufficient.
#'
#' A non-finite value on either side returns \code{FALSE}, so a diverged or
#' \code{NA} log-likelihood never terminates the EM loop by appearing
#' converged.
#'
#' @param old_val Criterion value at the previous iteration.
#' @param new_val Criterion value at the current iteration.
#' @param abs_tol Absolute-change tolerance (default \code{1e-6}).
#' @param rel_tol Relative-change tolerance, taken against \code{|old_val|}
#'   (default \code{1e-6}).
#' @return A single logical: \code{TRUE} if either criterion is met and both
#'   values are finite, otherwise \code{FALSE}.
#' @noRd
check_convergence <- function(old_val, new_val, abs_tol = 1e-6, rel_tol = 1e-6) {
  # A non-finite value must never read as converged: EM would stop on a
  # diverged log-likelihood and report the last valid parameters as the MLE.
  if (!is.finite(old_val) || !is.finite(new_val)) {
    return(FALSE)
  }

  abs_diff <- abs(new_val - old_val)
  rel_diff <- abs_diff / (abs(old_val) + .Machine$double.eps)

  abs_diff < abs_tol || rel_diff < rel_tol
}

#' Normalise a vector to sum to one, tolerating degenerate input
#'
#' Divides by the sum after flooring entries at \code{min_val}, which prevents
#' an exactly-zero posterior weight from removing an observation from a
#' cluster update irrecoverably. When the input contains a non-finite value, or
#' sums to something non-finite or non-positive, the uniform distribution
#' \eqn{1/n} is returned instead: in the E-step this expresses "no information
#' about this observation's cluster" rather than propagating \code{NaN} into
#' the M-step.
#'
#' @param x A numeric vector of non-negative weights. A zero-length vector is
#'   returned unchanged.
#' @param min_val Floor applied to each entry before normalising (default
#'   \code{.Machine$double.eps}).
#' @return A numeric vector the same length as \code{x} summing to one, or
#'   \code{x} unchanged if it has length zero.
#' @noRd
safe_normalize <- function(x, min_val = .Machine$double.eps) {
  .lucid_safe_normalize(x, min_val = min_val)
}


# =============================================================================
# Internal `.lucid_*` implementations (merged in from the former
# 00_stability_preload.R). Each exported wrapper above is a thin pass-through
# to one of these.
# =============================================================================

#' Log-sum-exp implementation behind \code{safe_log_sum_exp()}
#'
#' Row-wise if \code{x} is a matrix, otherwise scalar. Handles the
#' all-\code{-Inf}/all-\code{Inf} edge cases explicitly since
#' \code{max(x) + log(sum(exp(x - max(x))))} would otherwise produce \code{NaN}.
#'
#' @param x A numeric vector, or a matrix (log-sum-exp taken over each row).
#' @return A scalar if \code{x} is a vector, or a vector of one value per row
#'   if \code{x} is a matrix.
#' @noRd
.lucid_safe_log_sum_exp <- function(x) {
  if (is.null(dim(x))) {
    max_x <- max(x)
    if (is.infinite(max_x) && max_x < 0) return(-Inf)
    if (is.infinite(max_x) && max_x > 0) return(Inf)
    return(max_x + log(sum(exp(x - max_x))))
  }

  max_x <- apply(x, 1, max)
  out <- max_x + log(rowSums(exp(sweep(x, 1, max_x))))
  out[is.infinite(max_x) & max_x < 0] <- -Inf
  out[is.infinite(max_x) & max_x > 0] <- Inf
  out
}

#' Normalization implementation behind \code{safe_normalize()}
#'
#' @param x A numeric vector of non-negative weights.
#' @param min_val Floor applied to each entry before normalising.
#' @return \code{x} rescaled to sum to one, or the uniform distribution if
#'   \code{x} is non-finite, sums to a non-finite or non-positive value, or
#'   has length zero (unchanged in the last case).
#' @noRd
.lucid_safe_normalize <- function(x, min_val = .Machine$double.eps) {
  if (length(x) == 0) {
    return(x)
  }
  if (any(!is.finite(x))) {
    return(rep(1 / length(x), length(x)))
  }
  x <- pmax(x, min_val)
  s <- sum(x)
  if (!is.finite(s) || s <= 0) {
    return(rep(1 / length(x), length(x)))
  }
  x / s
}

#' Covariance-conditioning implementation behind \code{check_and_stabilize_sigma()}
#'
#' @param sigma A square covariance matrix.
#' @param threshold Condition-number threshold above which a ridge is added.
#' @param epsilon Ridge magnitude and eigenvalue/variance floor.
#' @return A symmetric, finite, well-conditioned matrix, or \code{sigma}
#'   unchanged if it is not square.
#' @noRd
.lucid_check_and_stabilize_sigma <- function(sigma, threshold = 1e10, epsilon = 1e-6) {
  if (!is.matrix(sigma) || nrow(sigma) != ncol(sigma)) {
    return(sigma)
  }
  # Sanitize non-finite entries before anything downstream sees them.  A NA or
  # Inf in the covariance used to be caught only incidentally, by whichever
  # density routine happened to fail on it.  Fall back to the (usually finite)
  # diagonal, which keeps row-wise variation in the density, and replace any
  # remaining non-finite variance with a small positive value.
  if (!all(is.finite(sigma))) {
    d <- diag(sigma)
    d[!is.finite(d) | d <= 0] <- epsilon
    sigma <- diag(d, nrow(sigma))
  }
  sigma <- (sigma + t(sigma)) / 2
  cond_num <- try(kappa(sigma), silent = TRUE)
  if (inherits(cond_num, "try-error") || !is.finite(cond_num) || cond_num > threshold) {
    sigma <- sigma + diag(epsilon, nrow(sigma))
  }
  eig <- try(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  if (inherits(eig, "try-error") || any(!is.finite(eig)) || min(eig) < epsilon) {
    sigma <- sigma + diag(epsilon, nrow(sigma))
  }
  sigma
}

#' Matrix-inversion implementation behind \code{safe_solve()}
#'
#' @param mat A square matrix to invert.
#' @param threshold Condition-number threshold passed to
#'   \code{.lucid_check_and_stabilize_sigma()} if a first \code{solve()} fails.
#' @param epsilon Ridge/floor magnitude passed through to the same fallback.
#' @return The inverse (or a pseudo-inverse, as a last resort) of \code{mat},
#'   or \code{NULL} if \code{mat} is not square or every fallback fails.
#' @noRd
.lucid_safe_solve <- function(mat, threshold = 1e10, epsilon = 1e-6) {
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    return(NULL)
  }

  result <- try(solve(mat), silent = TRUE)
  if (!inherits(result, "try-error")) {
    return(result)
  }

  stabilized <- .lucid_check_and_stabilize_sigma(mat, threshold = threshold, epsilon = epsilon)
  result <- try(solve(stabilized), silent = TRUE)
  if (!inherits(result, "try-error")) {
    return(result)
  }

  sv <- try(svd(stabilized), silent = TRUE)
  if (inherits(sv, "try-error")) {
    return(NULL)
  }
  d <- sv$d
  if (length(d) == 0 || all(d <= epsilon)) {
    return(NULL)
  }
  d_inv <- ifelse(d > epsilon, 1 / d, 0)
  sv$v %*% (diag(d_inv, nrow = length(d_inv)) %*% t(sv$u))
}


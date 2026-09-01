# Utility functions for numerical stability

# safe_log_sum_exp() is defined and documented in `stability_utils.R`.

#' Safe log operation with minimum value
#' 
#' @param x Numeric value
#' @param min_val Minimum value before taking log (default: 1e-300)
#' @return Safe log value
#' @noRd
safe_log <- function(x, min_val = 1e-300) {
  log(pmax(x, min_val))
}

# safe_solve() is defined and documented in `stability_utils.R`.

#' Log multivariate-normal density evaluated on observed coordinates only
#'
#' For each row of \code{Z}, evaluates the Gaussian log-density using only the
#' columns that are observed, i.e. the marginal density implied by the model.
#' A row with no observed value contributes 0, which is exactly the list-wise
#' partition of vbae123 Eq 11 (the Z term drops out of the likelihood).
#'
#' Rows are grouped by their missingness pattern so each distinct pattern
#' requires a single Cholesky factorization rather than one per row.
#'
#' @param Z numeric matrix, possibly containing NA
#' @param mean numeric vector of length ncol(Z)
#' @param sigma covariance matrix, ncol(Z) by ncol(Z)
#' @return numeric vector of length nrow(Z)
#' @noRd
log_dmvnorm_observed <- function(Z, mean, sigma) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  out <- numeric(n)
  if (n == 0L) return(out)

  # Guard the parameters explicitly.  Previously an NA in `mean` or `sigma` was
  # caught only incidentally, by mclust::dmvnorm raising an error; evaluating
  # the density directly would otherwise return NA silently.
  if (!all(is.finite(mean))) {
    stop("non-finite values in the cluster mean passed to the Gaussian density")
  }
  if (!all(is.finite(sigma))) {
    stop("non-finite values in the covariance matrix passed to the Gaussian density")
  }

  obs_mat <- !is.na(Z)
  # group rows sharing an identical observed-coordinate pattern
  keys <- apply(obs_mat, 1L, function(v) paste0(as.integer(v), collapse = ""))
  for (key in unique(keys)) {
    rows <- which(keys == key)
    cols <- which(obs_mat[rows[1L], ])
    if (length(cols) == 0L) {
      out[rows] <- 0
      next
    }
    s_sub <- sigma[cols, cols, drop = FALSE]
    s_sub <- check_and_stabilize_sigma(s_sub)
    ch <- try(chol(s_sub), silent = TRUE)
    if (inherits(ch, "try-error")) {
      # diagonal fallback: still a proper marginal, just an independent one
      sd_sub <- sqrt(pmax(diag(s_sub), .Machine$double.eps))
      out[rows] <- rowSums(stats::dnorm(
        Z[rows, cols, drop = FALSE],
        mean = matrix(mean[cols], nrow = length(rows), ncol = length(cols), byrow = TRUE),
        sd = matrix(sd_sub, nrow = length(rows), ncol = length(cols), byrow = TRUE),
        log = TRUE
      ))
      next
    }
    d <- length(cols)
    centered <- sweep(Z[rows, cols, drop = FALSE], 2L, mean[cols], "-")
    q <- backsolve(ch, t(centered), transpose = TRUE)
    quad <- colSums(q^2)
    out[rows] <- -0.5 * d * log(2 * pi) - sum(log(diag(ch))) - 0.5 * quad
  }
  out
}

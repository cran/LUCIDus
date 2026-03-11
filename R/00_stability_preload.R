# Preload core numerical helpers so files that alias them at top-level
# (e.g., lse <- safe_log_sum_exp) can be sourced regardless of collation order.

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

.lucid_check_and_stabilize_sigma <- function(sigma, threshold = 1e10, epsilon = 1e-6) {
  if (!is.matrix(sigma) || nrow(sigma) != ncol(sigma)) {
    return(sigma)
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

safe_log_sum_exp <- function(x) {
  .lucid_safe_log_sum_exp(x)
}

safe_normalize <- function(x, min_val = .Machine$double.eps) {
  .lucid_safe_normalize(x, min_val = min_val)
}

check_and_stabilize_sigma <- function(sigma, threshold = 1e10, epsilon = 1e-6) {
  .lucid_check_and_stabilize_sigma(sigma, threshold = threshold, epsilon = epsilon)
}

safe_solve <- function(mat, threshold = 1e10, epsilon = 1e-6) {
  .lucid_safe_solve(mat, threshold = threshold, epsilon = epsilon)
}

# Family helpers: canonical external labels are "normal"/"binary".
# Keep backward compatibility with "gaussian"/"binomial".
normalize_family_label <- function(family) {
  f <- tolower(as.character(family[[1]]))
  if (f %in% c("normal", "gaussian")) return("normal")
  if (f %in% c("binary", "binomial")) return("binary")
  stop("Unsupported family: ", family, ". Use one of normal/binary (gaussian/binomial supported as aliases).")
}

is_normal_family <- function(family) {
  normalize_family_label(family) == "normal"
}

is_binary_family <- function(family) {
  normalize_family_label(family) == "binary"
}

to_parallel_family <- function(family) {
  if (is_normal_family(family)) "gaussian" else "binomial"
}

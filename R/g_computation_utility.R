#' G-computation cluster log-likelihood, early model
#'
#' Like \code{\link{Estep_early}}, but the exposure (\eqn{G \to X}) term
#' only -- g-computation predicts from the exposure path alone, ignoring
#' \code{Z}/\code{Y} even if supplied.
#'
#' @param beta,mu,sigma,gamma Current parameter estimates (only \code{beta}
#'   is used).
#' @param G,Z,Y Exposure, omics, and outcome data (\code{Z}/\code{Y}
#'   accepted for a consistent call signature but unused).
#' @param family.list Unused; accepted for a consistent call signature.
#' @param K,N Number of clusters and observations.
#' @param useY,ind.na Unused; accepted for a consistent call signature.
#' @param ... Unused.
#' @return An \code{N x K} matrix of exposure-only log-likelihoods.
#' @noRd
Estep_early_g <- function(beta,
                          mu,
                          sigma,
                          gamma = NULL,
                          G,
                          Z,
                          Y = NULL,
                          family.list,
                          K,
                          N,
                          useY,
                          ind.na, ...) {
  # initialize vectors for storing likelihood
  pXgG <- pZgX <- pYgX <- matrix(rep(0, N * K), nrow = N)
  
  # log-likelihood for G -> X
  xb <- cbind(rep(1, N), G) %*% t(beta)
  xb_lse <- apply(xb, 1, lse)
  pXgG <- xb - xb_lse
  
  
  
  vec <- pXgG 
  return (vec)
}


#' G-computation joint-state log-likelihood, parallel model
#'
#' Like \code{\link{Estep}}, but the exposure (\eqn{G \to X}) term only,
#' summed across layers for each joint cluster-state combination.
#'
#' @param G,Z,Y Exposure, omics, and outcome data (\code{Z}/\code{Y}
#'   accepted for a consistent call signature but unused).
#' @param Beta Per-layer exposure coefficients.
#' @param Mu,Sigma,family,useY,na_pattern Unused; accepted for a consistent
#'   call signature.
#' @param Delta Outcome-model object; only \code{Delta$K} (per-layer cluster
#'   counts) is used.
#' @return A \code{K[1] x ... x K[nOmics] x N} array of exposure-only
#'   log-likelihoods.
#' @noRd
Estep_g <- function(G, Z, Y, Beta, Mu, Sigma, Delta, family, useY, na_pattern) {
  N <- nrow(G)
  K <- Delta$K
  f1 <- lapply(seq_along(K), function(i) f_GtoX(G, Beta[[i]]))
  states <- expand.grid(lapply(K, seq_len), KEEP.OUT.ATTRS = FALSE)
  out <- matrix(0, nrow = prod(K), ncol = N)
  for (s in seq_len(nrow(states))) {
    for (a in seq_along(K)) out[s, ] <- out[s, ] + f1[[a]][, states[s, a]]
  }
  array(out, dim = c(K, N))
}




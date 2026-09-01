#' Outcome-family dispatch object for a normal (continuous) Y, early model
#'
#' Returns the four closures the early EM loop needs to fit a continuous
#' outcome: random initialization (\code{initial.gamma}), the E-step density
#' (\code{f.pYgX}), the M-step (\code{f.maxY}, via
#' \code{\link{fit_early_outcome}}), and cluster relabeling
#' (\code{f.switch}, via \code{\link{relabel_early_parameters}}).
#'
#' @param K Number of latent clusters.
#' @param ... Unused; present so \code{normal}/\code{binary} share a call
#'   signature at the dispatch site.
#' @return A list: \code{initial.gamma}, \code{f.pYgX}, \code{f.maxY},
#'   \code{f.switch}.
#' @noRd
normal <- function(K, ...) {
  initial.gamma <- function(K, dimCoY) {
    early_gamma_from_levels(
      cluster_effect = stats::runif(K, min = -1, max = 1),
      covariate = stats::runif(dimCoY, min = -1, max = 1),
      sigma = stats::runif(K),
      covariate_names = if (dimCoY) paste0("CoY", seq_len(dimCoY)) else character(0)
    )
  }

  f.pYgX <- function(Y, gamma, K, N, CoY, dimCoY, itr) {
    levels <- early_gamma_levels(gamma, K)
    eta <- matrix(rep(levels, each = N), nrow = N)
    covariate <- early_gamma_covariate(gamma, K)
    if (dimCoY != 0 && length(covariate)) {
      eta <- eta + as.numeric(CoY %*% covariate)
    }
    sigma <- gamma$sigma
    vapply(seq_len(K), function(i) {
      stats::dnorm(as.numeric(Y), mean = eta[, i], sd = sigma[i], log = TRUE)
    }, numeric(N))
  }

  f.maxY <- function(Y, r, CoY, K, CoYnames) {
    fit_early_outcome(Y, r, CoY, K, CoYnames, family = "normal")
  }

  f.switch <- function(beta, mu, sigma, gamma, K, Tr = NULL) {
    relabel_early_parameters(beta, mu, sigma, gamma, K)
  }

  list(initial.gamma = initial.gamma, f.pYgX = f.pYgX,
       f.maxY = f.maxY, f.switch = f.switch)
}

#' Outcome-family dispatch object for a binary Y, early model
#'
#' Same role as \code{\link{normal}}, for a binary outcome fit with a
#' logistic \eqn{X \to Y} model.
#'
#' @param K Number of latent clusters.
#' @param ... Unused; present so \code{normal}/\code{binary} share a call
#'   signature at the dispatch site.
#' @return A list: \code{n.par}, \code{initial.gamma}, \code{f.pYgX},
#'   \code{f.maxY}, \code{f.switch}.
#' @noRd
binary <- function(K, ...) {
  initial.gamma <- function(K, dimCoY) {
    early_gamma_from_levels(
      cluster_effect = stats::runif(K, min = -1, max = 1),
      covariate = stats::runif(dimCoY, min = -1, max = 1),
      covariate_names = if (dimCoY) paste0("CoY", seq_len(dimCoY)) else character(0)
    )
  }

  f.pYgX <- function(Y, gamma, K, N, CoY, dimCoY, itr) {
    levels <- early_gamma_levels(gamma, K)
    eta <- matrix(rep(levels, each = N), nrow = N)
    covariate <- early_gamma_covariate(gamma, K)
    if (dimCoY != 0 && length(covariate)) {
      eta <- eta + as.numeric(CoY %*% covariate)
    }
    probability <- stats::plogis(eta)
    vapply(seq_len(K), function(i) {
      stats::dbinom(as.numeric(Y), 1, probability[, i], log = TRUE)
    }, numeric(N))
  }

  f.maxY <- function(Y, r, CoY, K, CoYnames) {
    fit_early_outcome(Y, r, CoY, K, CoYnames, family = "binary")
  }

  f.switch <- function(beta, mu, sigma, gamma, K) {
    relabel_early_parameters(beta, mu, sigma, gamma, K)
  }

  list(n.par = K, initial.gamma = initial.gamma, f.pYgX = f.pYgX,
       f.maxY = f.maxY, f.switch = f.switch)
}

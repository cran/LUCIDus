# General-purpose input validation and outcome-family labeling helpers used
# across the fitting, tuning, bootstrap and prediction entry points. Neither
# concern is specific to numerical stability (see `stability_utils.R`) or to
# any one model type.

#' Reject a non-omics input that contains missing or non-finite values
#'
#' Only the omics data may contain missing values -- that is what the
#' incomplete-omics extension models. Exposures, outcome and covariates must
#' be complete: a missing value there is not handled anywhere downstream, and
#' the symptom is not an error but a fit that either runs for a very long
#' time or returns badly wrong estimates with nothing in the output to
#' indicate it.
#'
#' @param x The input to check (matrix, vector, or \code{NULL}).
#' @param name Human-readable name of \code{x}, used in the error message.
#' @return Invisibly \code{NULL} if \code{x} is \code{NULL} or fully valid;
#'   otherwise raises an error naming the problem.
#' @noRd
check_complete_input <- function(x, name) {
  if (is.null(x)) return(invisible(NULL))
  if (anyNA(x)) {
    n <- sum(is.na(x))
    stop(name, " contains ", n, " missing value", if (n > 1) "s" else "",
         ". Only the omics data Z may contain NA; ",
         name, " must be complete. Remove or impute these observations first.",
         call. = FALSE)
  }
  if (is.numeric(x) && any(is.infinite(x))) {
    stop(name, " contains non-finite values (Inf or -Inf).", call. = FALSE)
  }
  invisible(NULL)
}

#' Normalize an outcome-family label to "normal" or "binary"
#'
#' Canonical external labels are "normal"/"binary"; "gaussian"/"binomial"
#' are accepted as aliases for backward compatibility.
#'
#' @param family Character, length 1 (or a vector, of which only the first
#'   element is used): one of "normal", "binary", "gaussian", "binomial".
#' @return "normal" or "binary".
#' @noRd
normalize_family_label <- function(family) {
  f <- tolower(as.character(family[[1]]))
  if (f %in% c("normal", "gaussian")) return("normal")
  if (f %in% c("binary", "binomial")) return("binary")
  stop("Unsupported family: ", family, ". Use one of normal/binary (gaussian/binomial supported as aliases).")
}

#' Is this outcome family normal?
#'
#' @param family A family label accepted by \code{normalize_family_label}.
#' @return Logical.
#' @noRd
is_normal_family <- function(family) {
  normalize_family_label(family) == "normal"
}

#' Is this outcome family binary?
#'
#' @param family A family label accepted by \code{normalize_family_label}.
#' @return Logical.
#' @noRd
is_binary_family <- function(family) {
  normalize_family_label(family) == "binary"
}

#' Translate a family label to the name the parallel model's GLM code expects
#'
#' The parallel model's outcome fitting goes through \code{glm()}-style family
#' names ("gaussian"/"binomial") rather than the package's own "normal"/
#' "binary" labels.
#'
#' @param family A family label accepted by \code{normalize_family_label}.
#' @return "gaussian" or "binomial".
#' @noRd
to_parallel_family <- function(family) {
  if (is_normal_family(family)) "gaussian" else "binomial"
}

#' @title Wrapper for LUCID Model and Penalty Tuning
#' @description Fit a grid of LUCID models over candidate numbers of latent
#' clusters \code{K} and (optionally) L1 penalties \code{Rho_G},
#' \code{Rho_Z_Mu}, and \code{Rho_Z_Cov}. The input format for \code{K} differs
#' by \code{lucid_model}. For "early", use an integer vector (for example,
#' \code{2:4}). For "parallel", use a list of vectors/integers, one per layer
#' (for example, \code{list(2:3, 2:3, 2)}). For "serial", use a nested list as
#' required by the serial model.
#' @param G Exposures, a numeric vector, matrix, or data frame. Categorical
#' variables should be transformed into dummy variables.
#' @param Z Omics data. If "early", an N by M matrix. If "parallel", a list of
#' matrices (same N). If "serial", a list matching the serial model structure.
#' @param Y Outcome, a numeric vector. Binary outcomes should be coded as 0/1.
#' @param CoG Optional covariates for the G-to-X model.
#' @param CoY Optional covariates for the X-to-Y model.
#' @param family Outcome family: "normal" or "binary".
#' @param K Candidate latent-cluster values in model-specific format.
#' @param lucid_model LUCID model type: "early", "parallel", or "serial".
#' @param Rho_G Scalar or vector penalty for exposure coefficients in the
#' G-to-X model. \code{CoG} covariates are included unpenalized. Vector tuning
#' is supported for "early" and "parallel". For "serial", only scalar inputs
#' are supported.
#' @param Rho_Z_Mu Scalar or vector penalty for cluster-specific Z means. Vector
#' tuning is supported for "early" and "parallel". For "serial", only scalar
#' inputs are supported.
#' @param Rho_Z_Cov Scalar or vector penalty for cluster-specific Z covariance
#' matrices. Vector tuning is supported for "early" and "parallel". For
#' "serial", only scalar inputs are supported.
#' @param verbose_tune Logical; print tuning progress if \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{estimate_lucid}}.
#' @return A list holding the tuning table, every fitted candidate, and the
#' selected model. The element names differ by \code{lucid_model}:
#' \itemize{
#'   \item "early": \code{tune_list}, \code{res_model}, \code{best_model}
#'   \item "parallel" and "serial": \code{tune_K}, \code{model_list},
#'     \code{model_opt}
#' }
#'
#' The tuning table has one row per grid point, in the order the candidates were
#' fitted, and the fitted-model list is aligned with it by position. Its columns
#' are the grid coordinates -- \code{K} (one column per layer or stage for
#' "parallel" and "serial") together with \code{Rho_G}, \code{Rho_Z_Mu} and
#' \code{Rho_Z_Cov} -- followed by \code{BIC}.
#'
#' Selection is on \code{BIC}, minimised over the rows. For "early" and
#' "parallel" this is the penalized BIC of Eq 18: the full parameter count is
#' reduced by one for each variable deselected by the penalty, so a sparser fit
#' is not charged for coefficients it has driven to zero. Penalties are not
#' tuned for "serial" -- only scalar penalties are accepted there -- so a serial
#' grid varies \code{K} alone. A candidate whose EM algorithm failed
#' records \code{NA} and is skipped; if every candidate fails, an error is
#' raised rather than a model returned. Ties are broken by taking the first
#' minimising row.
#'
#' Note that the returned optimum is the penalized fit itself. It is
#' \code{\link{lucid}}, not this function, that refits the selected variables
#' without a penalty -- so estimates taken straight from \code{best_model} or
#' \code{model_opt} here are shrunk towards zero.
#' @export
#' @examples
#' \dontrun{
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y <- sim_data$Y_normal
#' tune_early <- tune_lucid(G = G, Z = Z, Y = Y, lucid_model = "early", K = 2:3)
#' tune_rho <- tune_lucid(
#'   G = G, Z = Z, Y = Y, lucid_model = "early", K = 2,
#'   Rho_G = c(0, 0.1), Rho_Z_Mu = c(0, 5), Rho_Z_Cov = c(0, 0.1)
#' )
#' }
tune_lucid <-
function (G, Z, Y, CoG = NULL, CoY = NULL, family = c("normal", 
    "binary"), K, lucid_model = c("early", "parallel", "serial"), 
    Rho_G = 0, Rho_Z_Mu = 0, Rho_Z_Cov = 0, verbose_tune = FALSE, 
    ...) 
{
    family <- normalize_family_label(family)
    # Resolve once, as estimate_lucid() does for the same reason: forwarding
    # the raw, unresolved 3-element default to tune_lucid_auxi() (whose own
    # default is only c("early", "parallel")) made match.arg() reject it
    # there, so tune_lucid() could not actually be called without naming
    # lucid_model explicitly, despite its signature advertising a default.
    lucid_model <- match.arg(lucid_model)
    if (lucid_model == "early" || lucid_model == "parallel") {
        tune_lucid_auxi(G = G, Z = Z, Y = Y, CoG = CoG, CoY = CoY,
            family = family, K = K, lucid_model = lucid_model,
            Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu, Rho_Z_Cov = Rho_Z_Cov,
            verbose_tune = verbose_tune, ...)
    }
    else if (lucid_model == "serial") {
        if (length(Rho_G) > 1 | length(Rho_Z_Mu) > 1 | length(Rho_Z_Cov) > 
            1) {
            stop("Tune LUCID in serial can't tune for multiple regualrities for now!")
        }
        if (all(!sapply(K, is.list))) {
            K_list = expand_list_grid(K)
            for (i in 1:length(K_list)) {
                K_list[[i]] <- K_list[[i]][-length(K_list[[i]])]
                all_not_lists <- all(!sapply(K_list[[i]], is.list))
                if (all_not_lists) {
                  K_list[[i]] = unlist(K_list[[i]])
                }
            }
        }
        else {
            K_list = expand_list_grid(K)
            for (i in 1:length(K_list)) {
                K_list[[i]] <- K_list[[i]][-length(K_list[[i]])]
                for (j in 1:length(K_list[[i]])) {
                  if (length(K_list[[i]][[j]]) > 1) {
                    names(K_list[[i]][[j]]) = NULL
                    K_list[[i]][[j]] <- as.list(K_list[[i]][[j]])
                  }
                }
            }
        }
        K_matrix = as.matrix(K_list)
        m <- length(K_list)
        bic <- rep(NA_real_, m)
        model_list <- vector(mode = "list", length = m)
        if (verbose_tune && m > 1) {
            cat("Tuning LUCID model \n \n")
        }
        for (i in 1:m) {
            # Wrapped in try(), matching the early/parallel branches above:
            # one candidate K failing to converge (e.g. too many clusters for
            # the sample size) used to crash the whole tuning run instead of
            # being skipped like it is for early/parallel.
            fit <- try(estimate_lucid(G = G, Z = Z, Y = Y,
                CoG = CoG, CoY = CoY, family = family, lucid_model = "serial",
                K = K_list[[i]], Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu,
                Rho_Z_Cov = Rho_Z_Cov, verbose = verbose_tune,
                ...), silent = TRUE)
            if ("try-error" %in% class(fit)) {
                model_list[[i]] <- fit
                next
            }
            model_list[[i]] <- fit
            bic[i] <- cal_bic_serial(fit)
        }
        if (all(is.na(bic))) {
            stop("LUCID model fails to converge given current tuning parameters. ",
                "Last candidate's error: ", .last_try_error_message(model_list))
        }
        min_bic <- min(bic, na.rm = TRUE)
        # [1] guards against ties: which() can return several indices, and
        # `[[` with a length > 1 index does recursive indexing, not selection.
        model_opt_index <- which(bic == min_bic)[1]
        K_matrix <- cbind(K_matrix, bic)
        return(list(tune_K = K_matrix, model_list = model_list,
            model_opt = model_list[[model_opt_index]]))
    }
}

#' Grid search over K and penalties, for "early" and "parallel"
#'
#' Fits every combination of \code{K} and the penalty grids, records each
#' candidate's BIC (\code{NA} if it fails to converge, via \code{try()}), and
#' returns the lowest-BIC fit. \code{tune_lucid()}'s serial branch implements
#' the equivalent logic itself rather than calling this, since serial's
#' \code{K} has a different (nested-list) structure.
#'
#' @param G,Z,Y,CoG,CoY Exposure, omics, outcome, and optional covariate data.
#' @param family "normal" or "binary".
#' @param K Candidate cluster counts: a vector (early) or a list of per-layer
#'   vectors (parallel).
#' @param lucid_model "early" or "parallel".
#' @param Rho_G,Rho_Z_Mu,Rho_Z_Cov Candidate penalty grids.
#' @param verbose_tune Whether to print tuning progress.
#' @param ... Passed to \code{\link{estimate_lucid}}.
#' @return A list: \code{tune_K}/\code{tune_list} (the grid with BIC per
#'   row), \code{model_list} (every candidate fit, or a \code{try-error}),
#'   and \code{best_model}/\code{model_opt} (the lowest-BIC fit).
#' @noRd
tune_lucid_auxi <-
function (G, Z, Y, CoG = NULL, CoY = NULL, family = "normal",
    K = 2:5, lucid_model = c("early", "parallel"), Rho_G = 0,
    Rho_Z_Mu = 0, Rho_Z_Cov = 0, verbose_tune = FALSE, ...)
{
    family <- normalize_family_label(family)
    if (match.arg(lucid_model) == "early") {
        # Validate the whole K grid up front: without this, a bad candidate
        # (e.g. a typo'd K = 1) only surfaced once its own grid row failed
        # inside estimate_lucid() and was silently recorded as NA below,
        # rather than as a clear error about the tuning grid itself.
        check_K(K)
        tune_list <- expand.grid(K, Rho_G, Rho_Z_Mu, Rho_Z_Cov)
        colnames(tune_list) <- c("K", "Rho_G", "Rho_Z_Mu", "Rho_Z_Cov")
        m <- nrow(tune_list)
        tune_list$BIC <- rep(0, m)
        if (verbose_tune && m > 1) {
            cat("Tuning LUCID model \n \n")
        }
        res_model <- vector(mode = "list", length = m)
        for (i in 1:m) {
            fit <- try(estimate_lucid(G = G, Z = Z, Y = Y, CoG = CoG, 
                CoY = CoY, family = family, lucid_model = "early", 
                K = tune_list[i, 1], Rho_G = tune_list[i, 2], 
                Rho_Z_Mu = tune_list[i, 3], Rho_Z_Cov = tune_list[i, 
                  4], verbose = verbose_tune, ...), silent = TRUE)
            if ("try-error" %in% class(fit)) {
                tune_list[i, 5] <- NA
            }
            else {
                tune_list[i, 5] <- summary(fit)$BIC
            }
            res_model[[i]] <- fit
        }
        if (all(is.na(tune_list[, 5]))) {
            stop("LUCID model fails to converge given current tuning parameters. ",
                "Last candidate's error: ", .last_try_error_message(res_model))
        }
        x <- min(tune_list[, 5], na.rm = TRUE)
        # [1] guards against ties: which() can return several indices, and
        # `[[` with a length > 1 index does recursive indexing, not selection.
        best_model <- res_model[[which(tune_list[, 5] == x)[1]]]
        return(list(best_model = best_model, tune_list = tune_list, 
            res_model = res_model))
    }
    else if (match.arg(lucid_model) == "parallel") {
        K_list = K
        K_matrix <- as.matrix(expand.grid(K_list))
        if (min(K_matrix) < 2) {
            stop("minimum K should be 2")
        }
        rho_grid <- expand.grid(Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu,
            Rho_Z_Cov = Rho_Z_Cov)
        n_k <- nrow(K_matrix)
        n_rho <- nrow(rho_grid)
        n_total <- n_k * n_rho
        if (verbose_tune && n_total > 1) {
            cat("Tuning LUCID model \n \n")
        }
        bic <- rep(NA_real_, n_total)
        model_list <- vector(mode = "list", length = n_total)
        tune_rows <- vector(mode = "list", length = n_total)
        idx <- 0
        for (i in 1:n_k) {
            for (j in 1:n_rho) {
                idx <- idx + 1
                fit <- try(estimate_lucid(G = G, Z = Z, Y = Y, 
                  CoG = CoG, CoY = CoY, family = family, lucid_model = "parallel", 
                  K = as.numeric(K_matrix[i, ]), Rho_G = rho_grid$Rho_G[j], 
                  Rho_Z_Mu = rho_grid$Rho_Z_Mu[j], Rho_Z_Cov = rho_grid$Rho_Z_Cov[j], 
                  verbose = verbose_tune, ...), silent = TRUE)
                if (!("try-error" %in% class(fit))) {
                  bic[idx] <- cal_bic_parallel(fit)
                }
                model_list[[idx]] <- fit
                tune_rows[[idx]] <- c(as.numeric(K_matrix[i, ]), 
                  rho_grid$Rho_G[j], rho_grid$Rho_Z_Mu[j], rho_grid$Rho_Z_Cov[j], 
                  bic[idx])
            }
        }
        tune_K <- as.data.frame(do.call(rbind, tune_rows))
        k_names <- paste0("K", seq_len(ncol(K_matrix)))
        colnames(tune_K) <- c(k_names, "Rho_G", "Rho_Z_Mu", 
            "Rho_Z_Cov", "BIC")
        if (all(is.na(tune_K$BIC))) {
            stop("LUCID model fails to converge given current tuning parameters. ",
                "Last candidate's error: ", .last_try_error_message(model_list))
        }
        min_bic <- min(tune_K$BIC, na.rm = TRUE)
        model_opt_index <- which(tune_K$BIC == min_bic)[1]
        return(list(tune_K = tune_K, model_list = model_list, 
            model_opt = model_list[[model_opt_index]]))
    }
}

#' Pull the underlying error message out of the last try-error in a list
#'
#' `tune_lucid_auxi()`/`tune_lucid()`'s grid loops record a \code{try-error}
#' object for every candidate that fails, then swallow it into a generic
#' "fails to converge" message once every candidate has failed -- which used
#' to make a genuine EM non-convergence indistinguishable from an unrelated
#' error (e.g. a missing optional dependency) raised inside
#' \code{estimate_lucid()}. This surfaces the actual message from the last
#' failing candidate instead.
#'
#' @param model_list A list of fitted models and/or \code{try-error} objects.
#' @return The condition message of the last \code{try-error} in
#'   \code{model_list}, or \code{"(no error detail available)"} if none is
#'   found.
#' @noRd
.last_try_error_message <- function(model_list) {
    is_error <- vapply(model_list, function(x) inherits(x, "try-error"), logical(1))
    if (!any(is_error)) {
        return("(no error detail available)")
    }
    last_error <- model_list[is_error][[sum(is_error)]]
    msg <- conditionMessage(attr(last_error, "condition"))
    trimws(msg)
}

#' Cartesian product of nested tuning-grid values, preserving list structure
#'
#' Expands a list of candidate-value vectors (or nested lists of them, for a
#' serial model's per-stage grid) into one list per combination, recursing
#' from the last element inward.
#'
#' @param lst A list of candidate-value vectors or nested lists.
#' @return A list of combinations, each itself a list with one element per
#'   position in \code{lst}.
#' @noRd
expand_list_grid <-
function (lst)
{
    if (length(lst) == 0) {
        return(list(list(list())))
    }
    first_element <- lst[[1]]
    rest_of_list <- lst[-1]
    expanded_rest <- expand_list_grid(rest_of_list)
    grid <- expand.grid(first_element)
    result <- list()
    for (i in 1:nrow(grid)) {
        for (j in 1:length(expanded_rest)) {
            combined_element <- c(list(unlist(grid[i, ])), expanded_rest[[j]])
            result <- c(result, list(combined_element))
        }
    }
    return(result)
}

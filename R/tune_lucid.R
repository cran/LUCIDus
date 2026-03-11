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
#' @return A list containing the tuning table, fitted models, and the selected
#' optimal model. Returned element names differ slightly by model:
#' \itemize{
#'   \item "early": \code{best_model}, \code{tune_list}, \code{res_model}
#'   \item "parallel"/"serial": \code{model_opt}, \code{tune_K}, \code{model_list}
#' }
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
    if (match.arg(lucid_model) == "early" | match.arg(lucid_model) == 
        "parallel") {
        tune_lucid_auxi(G = G, Z = Z, Y = Y, CoG = CoG, CoY = CoY, 
            family = family, K = K, lucid_model = lucid_model, 
            Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu, Rho_Z_Cov = Rho_Z_Cov, 
            verbose_tune = verbose_tune, ...)
    }
    else if (match.arg(lucid_model) == "serial") {
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
        bic <- rep(0, length(K_list))
        model_list <- vector(mode = "list", length = length(K_list))
        for (i in 1:length(K_list)) {
            model_list[[i]] <- estimate_lucid(G = G, Z = Z, Y = Y, 
                CoG = CoG, CoY = CoY, family = family, lucid_model = "serial", 
                K = K_list[[i]], Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu, 
                Rho_Z_Cov = Rho_Z_Cov, verbose = verbose_tune, 
                ...)
            bic[i] <- cal_bic_serial(model_list[[i]])
        }
        model_opt_index <- which(bic == min(bic))
        K_matrix <- cbind(K_matrix, bic)
        return(list(tune_K = K_matrix, model_list = model_list, 
            model_opt = model_list[[model_opt_index]]))
    }
}

tune_lucid_auxi <-
function (G, Z, Y, CoG = NULL, CoY = NULL, family = "normal", 
    K = 2:5, lucid_model = c("early", "parallel"), Rho_G = 0, 
    Rho_Z_Mu = 0, Rho_Z_Cov = 0, verbose_tune = FALSE, ...) 
{
    family <- normalize_family_label(family)
    if (match.arg(lucid_model) == "early") {
        tune_list <- expand.grid(K, Rho_G, Rho_Z_Mu, Rho_Z_Cov)
        colnames(tune_list) <- c("K", "Rho_G", "Rho_Z_Mu", "Rho_Z_Cov")
        m <- nrow(tune_list)
        tune_list$BIC <- rep(0, m)
        if (m > 1) {
            cat("Tuning LUCID model \n \n")
        }
        res_model <- vector(mode = "list", length = m)
        for (i in 1:m) {
            fit <- try(estimate_lucid(G = G, Z = Z, Y = Y, CoG = CoG, 
                CoY = CoY, family = family, lucid_model = "early", 
                K = tune_list[i, 1], Rho_G = tune_list[i, 2], 
                Rho_Z_Mu = tune_list[i, 3], Rho_Z_Cov = tune_list[i, 
                  4], verbose = verbose_tune, ...))
            if ("try-error" %in% class(fit)) {
                tune_list[i, 5] <- NA
            }
            else {
                tune_list[i, 5] <- summary(fit)$BIC
            }
            res_model[[i]] <- fit
        }
        if (all(is.na(tune_list[, 5]))) {
            stop("LUCID model fails to converge given current tuning parameters")
        }
        x <- min(tune_list[, 5], na.rm = TRUE)
        best_model <- res_model[[which(tune_list[, 5] == x)]]
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
            stop("LUCID model fails to converge given current tuning parameters")
        }
        min_bic <- min(tune_K$BIC, na.rm = TRUE)
        model_opt_index <- which(tune_K$BIC == min_bic)[1]
        return(list(tune_K = tune_K, model_list = model_list, 
            model_opt = model_list[[model_opt_index]]))
    }
}

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

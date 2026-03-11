#' @title Fit LUCID models with one or multiple omics layers
#' @description EM algorithm to estimate LUCID with one or multiple omics layers
#' @param lucid_model Specifying LUCID model, "early" for early integration, "parallel" for lucid in parallel,
#' "serial" for lucid in serial
#' @param G an N by P matrix representing exposures
#' @param Z Omics data, if "early", an N by M matrix; If "parallel", a list, each element i is a matrix with N rows and P_i features;
#' If "serial", a list, each element i is a matrix with N rows and p_i features or a list with two or more matrices with N rows and a certain number of features
#' @param Y a length N vector
#' @param CoG an N by V matrix representing covariates to be adjusted for G -> X
#' @param CoY an N by K matrix representing covariates to be adjusted for X -> Y
#' @param K Number of latent clusters. If "early", an integer greater or equal to 2; If "parallel", an integer vector, same length as Z, with each element being an integer greater or equal to 2;
#' If "serial", a list, each element is either an integer like that for "early" or an list of integers like that for "parallel", same length as Z
#' @param init_omic.data.model a vector of strings specifies the geometric model of omics
#' data. If NULL, See more in ?mclust::mclustModelNames
#' @param useY logical, if TRUE, EM algorithm fits a supervised LUCID; otherwise
#' unsupervised LUCID.
#' @param tol stopping criterion for the EM algorithm
#' @param max_itr Maximum iterations of the EM algorithm. If the EM algorithm iterates
#' more than max_itr without converging, the EM algorithm is forced to stop.
#' @param max_tot.itr Max number of total iterations for \code{estimate_lucid} function.
#' \code{estimate_lucid} may conduct EM algorithm for multiple times if the algorithm
#' fails to converge.
#' @param Rho_G A scalar. This parameter is the LASSO penalty to regularize
#' exposure coefficients in the G-to-X model. \code{CoG} adjustment covariates
#' are included unpenalized. If user wants to tune the penalty, use the wrapper
#' function \code{lucid}. Penalty tuning is supported for "early" and
#' "parallel". For "serial", only scalar penalty inputs are supported.
#' @param Rho_Z_Mu A scalar. This parameter is the LASSO penalty to
#' regularize cluster-specific means for omics data (Z). If user wants to tune the
#' penalty, use the wrapper function \code{lucid}. Penalty tuning is supported
#' for "early" and "parallel". For "serial", only scalar penalty inputs are
#' supported.
#' @param Rho_Z_Cov A scalar. This parameter is the graphical LASSO
#' penalty to estimate sparse cluster-specific variance-covariance matrices for omics
#' data (Z). If user wants to tune the penalty, use the wrapper function \code{lucid}.
#' Penalty tuning is supported for "early" and "parallel". For "serial", only
#' scalar penalty inputs are supported.
#' @param family The distribution of the outcome
#' @param seed Random seed to initialize the EM algorithm
#' @param init_impute Method to initialize the imputation of missing values in
#' LUCID. \code{mix} will use \code{mclust:imputeData} to implement EM Algorithm
#' for Unrestricted General Location Model by the mix package to impute the missing values in omics
#' data; \code{lod} will initialize the imputation via replacing missing values by
#' LOD / sqrt(2). LOD is determined by the minimum of each variable in omics data.
#' @param init_par For "early", an interface to initialize EM algorithm, if mclust,
#' initiate the parameters using the \code{mclust} package, if random, initiate the parameters
#' by drawing from a uniform distribution;
#' For "parallel", mclust is the default for quick convergence;
#' For "serial", each sub-model follows the above depending on it is a "early" or "parallel"
#' @param verbose Logging level for fitting progress. If \code{FALSE}, concise
#' start/finish status lines are printed. If \code{TRUE}, detailed iteration-level
#' traces (including log-likelihood updates) are printed.
#'
#' @import mclust
#' @import stats
#' @import utils
#' @import glasso
#' @import glmnet
#' @return A list contains the object below:
#' 1. res_Beta: estimation for G->X associations
#' 2. res_Mu: estimation for the mu of the X->Z associations
#' 3. res_Sigma: estimation for the sigma of the X->Z associations
#' 4. res_Gamma: estimation for X->Y associations
#' 5. inclusion.p: inclusion probability of cluster assignment for each observation
#' 6. K: number of latent clusters for "early"/list of numbers of latent clusters for "parallel" and "serial"
#' 7. var.names: names for the G, Z, Y variables
#' 8. init_omic.data.model: pre-specified geometric model of multi-omics data
#' 9. likelihood: converged LUCID model log likelihood
#' 10. family: the distribution of the outcome
#' 11. select: for "early" and "parallel", feature-selection indicators.
#' For "parallel", \code{select$selectG} is the exposure-wise union across layers
#' (selected in at least one layer) and \code{select$selectG_layer} stores
#' per-layer exposure selection.
#' 12. useY: whether this LUCID model is supervised
#' 13. Z: multi-omics data
#' 14. init_impute: pre-specified imputation method
#' 15. init_par: pre-specified parameter initialization method
#' 16. Rho: for "early" and "parallel", pre-specified regularity tuning parameters
#' 17. N: number of observations
#' 18. submodel: for LUCID in serial only, storing all the submodels
#' 19. em_control: EM stopping controls used for fitting (\code{tol},
#' \code{max_itr}, \code{max_tot.itr}); used by bootstrap refits
#' @examples
#' i <- 1008
#' set.seed(i)
#' G <- matrix(rnorm(500), nrow = 100)
#' Z1 <- matrix(rnorm(1000),nrow = 100)
#' Z2 <- matrix(rnorm(1000), nrow = 100)
#' Z3 <- matrix(rnorm(1000), nrow = 100)
#' Z4 <- matrix(rnorm(1000), nrow = 100)
#' Z5 <- matrix(rnorm(1000), nrow = 100)
#' Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5)
#' Y <- rnorm(100)
#' CoY <- matrix(rnorm(200), nrow = 100)
#' CoG <- matrix(rnorm(200), nrow = 100)
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2,2,2,2,2),
#' lucid_model = "serial",
#' family = "normal",
#' seed = i,
#' CoG = CoG, CoY = CoY,
#' useY = TRUE)
#' @export
#'
estimate_lucid <- function(lucid_model = c("early", "parallel","serial"),
                      G, Z, Y, CoG = NULL, CoY = NULL, K,
                      init_omic.data.model = "EEV",
                      useY = TRUE,
                      tol = 1e-3,
                      max_itr = 1e3,
                      max_tot.itr = 1e4,
                      Rho_G = 0,
                      Rho_Z_Mu = 0,
                      Rho_Z_Cov = 0,
                      family = c("normal", "binary"),
                      seed = 123,
                      init_impute = c("mix", "lod"),
                      init_par = c("mclust", "random"),
                      verbose = FALSE) {
  family <- normalize_family_label(family)
  if (match.arg(lucid_model) == "early" | match.arg(lucid_model) == "parallel"){
    # ========================== Early Integration ==========================
    # ========================== LUCID IN PARALLEL ==========================
    results <- est_lucid(lucid_model = lucid_model,
                         G = G,
                         Z = Z,
                         Y = Y,
                         CoG = CoG,
                         CoY = CoY, K = K,
                         init_omic.data.model = init_omic.data.model,
                         useY = useY,
                         tol = tol,
                         max_itr = max_itr,
                         max_tot.itr = max_tot.itr,
                         Rho_G = Rho_G,
                         Rho_Z_Mu = Rho_Z_Mu,
                         Rho_Z_Cov = Rho_Z_Cov,
                         family = family,
                         seed = seed,
                         init_impute = init_impute,
                         init_par = init_par,
                         verbose = verbose)
    return(results)
  }else{
    # ========================== LUCID IN Serial ==========================
    if (!is.list(K)) {
      if (is.numeric(K) && length(K) > 0) {
        # Backward-compatible support: serial with only early stages may pass K as numeric vector.
        K <- as.list(as.numeric(K))
      } else {
        stop("For LUCID in Serial, K should be a non-empty list!")
      }
    }
    if (length(K) == 0) {
      stop("For LUCID in Serial, K should be a non-empty list!")
    }
    if (is.null(Z)) {
      stop("Input data 'Z' is missing")
    }
    if (!is.list(Z)) {
      stop("Input data 'Z' should be a list for LUCID in Serial!")
    }
    if (length(Z) != length(K)) {
      stop("Z and K should be two lists of the same length for LUCID in Serial!")
    }

    normalize_serial_block <- function(z_block, k_block, block_id = "root") {
      if (is.list(k_block)) {
        if (!is.list(z_block)) {
          stop(paste0("For LUCID in Serial, input data 'Z' must match 'K' structure. Nested K at ",
                      block_id, " requires nested Z list."))
        }
        if (length(z_block) != length(k_block)) {
          stop(paste0("For LUCID in Serial, nested Z length must equal nested K length at ",
                      block_id, "."))
        }
        out <- vector("list", length(z_block))
        for (j in seq_along(z_block)) {
          out[[j]] <- normalize_serial_block(z_block[[j]], k_block[[j]],
                                             block_id = paste0(block_id, ".", j))
        }
        names(out) <- names(z_block)
        return(out)
      }

      if (!is.numeric(k_block) || length(k_block) != 1 || is.na(k_block) || k_block < 2) {
        stop(paste0("For estimate_lucid(..., lucid_model = 'serial'), each non-list K entry ",
                    "must be a single integer >= 2. Invalid entry at ", block_id, "."))
      }

      z_mat <- as.matrix(z_block)
      if (!is.numeric(z_mat)) {
        stop(paste0("For LUCID in Serial, Z block at ", block_id, " must be numeric."))
      }
      z_mat
    }

    for (i in seq_along(K)) {
      Z[[i]] <- normalize_serial_block(Z[[i]], K[[i]], block_id = paste0("stage", i))
    }

    extract_nonref_pip <- function(model_obj, n_obs) {
      if (inherits(model_obj, "early_lucid")) {
        p <- as.matrix(model_obj$inclusion.p)
        if (ncol(p) <= 1) {
          return(matrix(numeric(0), nrow = n_obs))
        }
        return(p[, -1, drop = FALSE])
      }
      if (inherits(model_obj, "lucid_parallel")) {
        p_list <- model_obj$inclusion.p
        p_nonref <- lapply(p_list, function(p) {
          p <- as.matrix(p)
          if (ncol(p) <= 1) {
            matrix(numeric(0), nrow = nrow(p))
          } else {
            p[, -1, drop = FALSE]
          }
        })
        if (length(p_nonref) == 0) {
          return(matrix(numeric(0), nrow = n_obs))
        }
        return(do.call(cbind, p_nonref))
      }
      stop("Unsupported submodel class in serial pipeline.")
    }

    fit_serial_stage <- function(stage_idx,
                                 stage_model,
                                 G_stage,
                                 Z_stage,
                                 Y_stage,
                                 CoG_stage,
                                 CoY_stage,
                                 K_stage,
                                 useY_stage,
                                 family_stage) {
      G_stage <- as.matrix(G_stage)
      if (!is.numeric(G_stage)) {
        stop("Serial stage input G must be numeric.")
      }

      rho_g_stage <- Rho_G
      if (rho_g_stage > 0 && ncol(G_stage) < 2) {
        rho_g_stage <- 0
        if (verbose) {
          cat(sprintf("Sub Model %d: Rho_G reset to 0 because stage G has fewer than 2 variables.\n",
                      stage_idx))
        }
      }

      if (isTRUE(verbose)) {
        return(est_lucid(
          lucid_model = stage_model,
          G = G_stage,
          Z = Z_stage,
          Y = Y_stage,
          CoG = CoG_stage,
          CoY = CoY_stage,
          K = K_stage,
          init_omic.data.model = init_omic.data.model,
          useY = useY_stage,
          tol = tol,
          max_itr = max_itr,
          max_tot.itr = max_tot.itr,
          Rho_G = rho_g_stage,
          Rho_Z_Mu = Rho_Z_Mu,
          Rho_Z_Cov = Rho_Z_Cov,
          family = family_stage,
          seed = seed + stage_idx * 1900,
          init_impute = init_impute,
          init_par = init_par,
          verbose = TRUE
        ))
      }

      stage_fit <- NULL
      invisible(capture.output(
        stage_fit <- est_lucid(
          lucid_model = stage_model,
          G = G_stage,
          Z = Z_stage,
          Y = Y_stage,
          CoG = CoG_stage,
          CoY = CoY_stage,
          K = K_stage,
          init_omic.data.model = init_omic.data.model,
          useY = useY_stage,
          tol = tol,
          max_itr = max_itr,
          max_tot.itr = max_tot.itr,
          Rho_G = rho_g_stage,
          Rho_Z_Mu = Rho_Z_Mu,
          Rho_Z_Cov = Rho_Z_Cov,
          family = family_stage,
          seed = seed + stage_idx * 1900,
          init_impute = init_impute,
          init_par = init_par,
          verbose = FALSE
        )
      ))
      stage_fit
    }

    n_stage <- length(K)
    post.p.list <- vector(mode = "list", length = n_stage)
    res.mu.list <- vector(mode = "list", length = n_stage)
    res.sigma.list <- vector(mode = "list", length = n_stage)
    res.delta.list <- vector(mode = "list", length = max(0, n_stage - 1))
    Znames <- vector(mode = "list", length = n_stage)
    submodel <- vector(mode = "list", length = n_stage)
    missing_by_stage <- vector(mode = "list", length = n_stage)
    has_penalty <- function(model_obj) {
      if (is.null(model_obj$Rho)) return(FALSE)
      any(unlist(model_obj$Rho[c("Rho_G", "Rho_Z_Mu", "Rho_Z_Cov")]) != 0)
    }
    stage_selection_msg <- function(model_obj) {
      if (inherits(model_obj, "early_lucid")) {
        g_sel <- sum(model_obj$select$selectG)
        g_tot <- length(model_obj$select$selectG)
        z_sel <- sum(model_obj$select$selectZ)
        z_tot <- length(model_obj$select$selectZ)
        return(sprintf("Selected G: %d/%d; Selected Z: %d/%d.", g_sel, g_tot, z_sel, z_tot))
      }
      if (inherits(model_obj, "lucid_parallel")) {
        g_sel <- sum(model_obj$select$selectG)
        g_tot <- length(model_obj$select$selectG)
        z_sel <- vapply(model_obj$select$selectZ, function(x) {
          if (is.null(dim(x))) sum(x) else sum(colSums(x) > 0)
        }, numeric(1))
        z_tot <- vapply(model_obj$select$selectZ, function(x) {
          if (is.null(dim(x))) length(x) else ncol(x)
        }, numeric(1))
        return(sprintf("Selected G: %d/%d; Selected Z by layer: %s.",
                       g_sel, g_tot, paste0(z_sel, "/", z_tot, collapse = ", ")))
      }
      ""
    }

    if (!isTRUE(verbose)) {
      cat(sprintf("Fitting LUCID serial model (%d stages)...\n", n_stage))
    }

    post.p <- NULL
    for (stage_idx in seq_len(n_stage)) {
      if (verbose) {
        cat("Fitting LUCID in Serial model",
            paste0("(", "Sub Model Number = ", stage_idx, ")"),
            "\n")
      }

      is_last <- (stage_idx == n_stage)
      stage_model <- if (is.list(K[[stage_idx]])) "parallel" else "early"
      stage_K <- if (is.list(K[[stage_idx]])) {
        as.numeric(unlist(K[[stage_idx]], use.names = FALSE))
      } else {
        as.numeric(K[[stage_idx]])
      }
      stage_Y <- if (is_last) Y else runif(nrow(G))
      stage_family <- if (is_last) family else "normal"
      stage_useY <- if (is_last) useY else FALSE
      stage_CoY <- if (is_last) CoY else NULL
      stage_CoG <- if (stage_idx == 1) CoG else NULL
      stage_G <- if (stage_idx == 1) G else post.p

      temp_model <- fit_serial_stage(
        stage_idx = stage_idx,
        stage_model = stage_model,
        G_stage = stage_G,
        Z_stage = Z[[stage_idx]],
        Y_stage = stage_Y,
        CoG_stage = stage_CoG,
        CoY_stage = stage_CoY,
        K_stage = stage_K,
        useY_stage = stage_useY,
        family_stage = stage_family
      )
      if (!isTRUE(verbose)) {
        if (has_penalty(temp_model)) {
          cat(sprintf("  Stage %d/%d (%s) finished. %s\n",
                      stage_idx, n_stage, stage_model, stage_selection_msg(temp_model)))
        } else {
          cat(sprintf("  Stage %d/%d (%s) finished.\n",
                      stage_idx, n_stage, stage_model))
        }
      }

      post.p.list[[stage_idx]] <- temp_model$inclusion.p
      res.mu.list[[stage_idx]] <- temp_model$res_Mu
      res.sigma.list[[stage_idx]] <- temp_model$res_Sigma
      Znames[[stage_idx]] <- temp_model$var.names$Znames
      submodel[[stage_idx]] <- temp_model
      missing_by_stage[[stage_idx]] <- temp_model$missing_summary

      if (stage_idx == 1) {
        res_Beta <- temp_model$res_Beta
        Gnames <- temp_model$var.names$Gnames
      } else {
        res.delta.list[[stage_idx - 1]] <- temp_model$res_Beta
      }

      if (is_last) {
        res_Gamma <- temp_model$res_Gamma
        Ynames <- temp_model$var.names$Ynames
      } else {
        post.p <- extract_nonref_pip(temp_model, n_obs = nrow(G))
        if (ncol(post.p) == 0) {
          stop(paste0("Sub Model ", stage_idx,
                      " produced no non-reference cluster probabilities to pass forward."))
        }
      }
    }

    if(verbose){
      cat("Success: LUCID in Serial Model is constructed!", "\n\n")
    } else {
      cat("Finished LUCID serial model.\n")
    }
    serial_missing_summary <- list(
      n_stages = n_stage,
      stage = missing_by_stage
    )
    results <- list(res_Beta = res_Beta,
                    res_Mu = res.mu.list,
                    res_Sigma = res.sigma.list,
                    res_Delta = res.delta.list,
                    res_Gamma = res_Gamma,
                    K = K,
                    N = nrow(G),
                    var.names =list(Gnames = Gnames,
                                    Znames = Znames,
                                    Ynames = Ynames),
                    init_omic.data.model =  init_omic.data.model,
                    #likelihood = loglik_update, *??? needs to discuss
                    inclusion.p = post.p.list,
                    family = family,
                    #select = list(selectG = selectG, selectZ = selectZ),
                    useY = useY,
                    Z = Z,
                    #z = Estep_r,
                    init_impute = init_impute,
                    init_par = init_par,
                    submodel = submodel,
                    missing_summary = serial_missing_summary,
                    Rho = list(Rho_G = Rho_G,
                               Rho_Z_Mu = Rho_Z_Mu,
                               Rho_Z_Cov = Rho_Z_Cov),
                    em_control = list(tol = tol, max_itr = max_itr,
                                      max_tot.itr = max_tot.itr)
    )
    class(results) <- c("lucid_serial")
    return(results)
  }
}

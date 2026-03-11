#' Summarize results of the early LUCID model
#'
#' @param object A LUCID model fitted by \code{\link{estimate_lucid}}
#' @param ... Additional argument \code{boot.se}, which can be an object
#' returned by \code{\link{boot_lucid}} to display bootstrap CIs in print output.
#' @return A list containing model information, fit statistics, feature-selection
#' summaries, detailed parameter estimates, missing-data summaries, and optional
#' bootstrap CI tables.
#' @export
#' @examples 
#' \donttest{
#' # use simulated data
#' G <- sim_data$G[1:300, , drop = FALSE]
#' Z <- sim_data$Z[1:300, , drop = FALSE]
#' Y_normal <- sim_data$Y_normal[1:300]
#'
#' # fit lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", family = "normal", K = 2,
#' seed = 1008)
#'
#' # conduct bootstrap resampling
#' boot1 <- suppressWarnings(
#'   boot_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", model = fit1, R = 5)
#' )
#'
#' # summarize lucid model
#' summary(fit1)
#'
#' # summarize lucid model with bootstrap CIs
#' summary(fit1, boot.se = boot1)
#' }
summary_lucid <- function(object, ...) {
  summary(object, auto_print = FALSE, ...)
}

#' @rdname summary_lucid
#' @method summary early_lucid
#' @export
summary.early_lucid <- function(object, ...) {
  args <- list(...)
  auto_print <- if (is.null(args$auto_print)) TRUE else isTRUE(args$auto_print)
  
  # Extract feature selection results
  selectG <- object$select$selectG
  selectZ <- object$select$selectZ
  
  # Count selected features
  nG <- if(is.null(selectG)) ncol(object$G) else sum(selectG)
  nZ <- if(is.null(selectZ)) ncol(object$Z) else sum(selectZ)
  
  K <- object$K
  gamma <- object$res_Gamma$beta
  
  # Count parameters
  if(is_normal_family(object$family)){
    nY <- length(object$res_Gamma$beta) + length(object$res_Gamma$sigma)
  }
  if(is_binary_family(object$family)){
    nY <- length(object$res_Gamma$beta)
  }
  
  npars <- (nG + 1) * (K - 1) + (nZ * K + nZ^2 * K) + nY
  BIC <- -2 * object$likelihood + npars * log(nrow(object$inclusion.p))
  
  # Prepare feature selection summary
  feature_summary <- list()
  
  # G features summary
  if(!is.null(selectG)) {
    # Ensure lengths match by taking the shorter length
    n_features <- min(length(object$var.names$Gnames), length(selectG))
    G_features <- data.frame(
      Feature = object$var.names$Gnames[1:n_features],
      Selected = selectG[1:n_features],
      row.names = NULL
    )
    feature_summary$G <- G_features
  }
  
  # Z features summary
  if(!is.null(selectZ)) {
    # Ensure lengths match
    n_features <- min(length(object$var.names$Znames), length(selectZ))
    Z_features <- data.frame(
      Feature = object$var.names$Znames[1:n_features],
      Selected = selectZ[1:n_features],
      row.names = NULL
    )
    feature_summary$Z <- Z_features
  }
  
  # Prepare regularization summary if available
  reg_summary <- NULL
  if(!is.null(object$Rho)) {
    reg_summary <- list(
      Rho_G = object$Rho$Rho_G,
      Rho_Z_Mu = object$Rho$Rho_Z_Mu,
      Rho_Z_Cov = object$Rho$Rho_Z_Cov
    )
  }
  
  results <- list(
    BIC = BIC,
    loglik = object$likelihood,
    model_info = list(
      family = object$family,
      K = K,
      n_observations = nrow(object$inclusion.p),
      n_features = list(
        G = nG,
        Z = nZ
      )
    ),
    feature_selection = feature_summary,
    regularization = reg_summary,
    model_fit = list(
      BIC = BIC,
      loglik = object$likelihood,
      n_parameters = npars
    ),
    parameters = list(
      beta = object$res_Beta[, c(TRUE, selectG)],
      mu = object$res_Mu[, selectZ],
      gamma = object$res_Gamma
    ),
    boot.se = args$boot.se,
    missing_data = object$missing_summary
  )
  
  class(results) <- "sumlucid_early"
  if (auto_print) {
    print(results)
  }
  invisible(results)
}






#' Summarize results of the parallel LUCID model
#'
#' @param object A LUCID model fitted by \code{\link{estimate_lucid}}
#' @param ... Additional argument \code{boot.se}, which can be an object
#' returned by \code{\link{boot_lucid}} to display bootstrap CIs in print output.
#'
#' @return A list containing model information, fit statistics, feature-selection
#' summaries (overall + by layer), detailed parameter estimates (by layer),
#' missing-data summaries (by layer), and optional bootstrap CI tables.
#' @export
summary.lucid_parallel <- function(object, ...) {
  args <- list(...)
  auto_print <- if (is.null(args$auto_print)) TRUE else isTRUE(args$auto_print)
  # Get dimensions and parameters
  K <- object$K
  nOmics <- length(K)
  
  # Extract feature selection results
  selectG <- object$select$selectG
  selectG_layer <- object$select$selectG_layer
  selectZ <- object$select$selectZ
  if(is.null(selectG_layer) && is.list(selectG)) {
    selectG_layer <- selectG
  }
  if(is.null(selectG) && !is.null(selectG_layer) && length(selectG_layer) > 0) {
    selectG <- Reduce("|", selectG_layer)
  }
  if(is.list(selectG)) {
    selectG <- Reduce("|", selectG)
  }
  
  # Count selected features
  nG <- if(is.null(selectG)) {
    length(object$var.names$Gnames)
  } else {
    sum(selectG)
  }
  nG_layer <- if(is.null(selectG_layer)) {
    rep(nG, nOmics)
  } else {
    sapply(selectG_layer, sum)
  }
  if(length(nG_layer) != nOmics) {
    nG_layer <- rep(nG, nOmics)
  }
  
  nZ <- if(is.null(selectZ)) {
    sapply(object$Z, ncol)
  } else {
    sapply(selectZ, function(x) {
      if (is.null(dim(x))) {
        sum(x)
      } else {
        sum(colSums(x) > 0)
      }
    })
  }
  
  # Count parameters
  if(to_parallel_family(object$family) == "gaussian") {
    nY <- length(object$res_Gamma$Gamma$mu) + length(object$res_Gamma$Gamma$sd)
  } else if(to_parallel_family(object$family) == "binomial") {
    nY <- length(object$res_Gamma$fit$coefficients)
  }
  
  # Calculate total number of parameters
  npars <- 0
  # G to X parameters (including intercepts)
  for(i in 1:nOmics) {
    npars <- npars + (nG_layer[i] + 1) * (K[i] - 1)
  }
  # X to Z parameters (means and covariances)
  for(i in 1:nOmics) {
    npars <- npars + (nZ[i] * K[i] + nZ[i] * nZ[i] * K[i])
  }
  # X to Y parameters
  npars <- npars + nY
  
  # Calculate BIC
  N <- nrow(object$inclusion.p[[1]])
  BIC <- -2 * object$likelihood + npars * log(N)
  
  # Prepare feature selection summary
  feature_summary <- list()
  
  # G features summary (overall: selected in any layer)
  if(!is.null(selectG)) {
    n_features <- min(length(object$var.names$Gnames), length(selectG))
    feature_summary$G <- data.frame(
      Feature = object$var.names$Gnames[1:n_features],
      Selected = selectG[1:n_features],
      row.names = NULL
    )
  }

  # G features summary by layer
  if(!is.null(selectG_layer)) {
    G_features_layer <- vector("list", nOmics)
    for(i in 1:nOmics) {
      selected_i <- selectG_layer[[i]]
      n_features <- min(length(object$var.names$Gnames), length(selected_i))
      G_features_layer[[i]] <- data.frame(
        Feature = object$var.names$Gnames[1:n_features],
        Selected = selected_i[1:n_features],
        row.names = NULL
      )
    }
    feature_summary$G_layer <- G_features_layer
  }
  
  # Z features summary per layer
  if(!is.null(selectZ)) {
    Z_features <- vector("list", nOmics)
    for(i in 1:nOmics) {
      selected_i <- selectZ[[i]]
      if (is.null(dim(selected_i))) {
        selected_counts <- as.integer(selected_i)
        selected_any <- as.logical(selected_i)
      } else {
        selected_counts <- colSums(selected_i)
        selected_any <- selected_counts > 0
      }
      Z_features[[i]] <- data.frame(
        Feature = object$var.names$Znames[[i]],
        Selected_in_clusters = selected_counts,
        Selected = selected_any,
        row.names = NULL
      )
    }
    feature_summary$Z <- Z_features
  }
  
  # Prepare regularization summary if available
  reg_summary <- NULL
  if(!is.null(object$Rho)) {
    reg_summary <- list(
      Rho_G = object$Rho$Rho_G,
      Rho_Z_Mu = object$Rho$Rho_Z_Mu,
      Rho_Z_Cov = object$Rho$Rho_Z_Cov
    )
  }
  
  # Combine results
  results <- list(
    model_info = list(
      family = object$family,
      K = K,
      n_observations = N,
      n_features = list(
        G = nG,
        Z = nZ
      )
    ),
    feature_selection = feature_summary,
    regularization = reg_summary,
    model_fit = list(
      BIC = BIC,
      loglik = object$likelihood,
      n_parameters = npars
    ),
    parameters = list(
      beta = object$res_Beta,
      mu = object$res_Mu,
      Gamma = object$res_Gamma
    ),
    boot.se = args$boot.se,
    missing_data = object$missing_summary
  )
  
  class(results) <- "sumlucid_parallel"
  if (auto_print) {
    print(results)
  }
  invisible(results)
}


#' Summarize results of the serial LUCID model
#'
#' @param object A LUCID model fitted by \code{\link{estimate_lucid}}
#' @param ... Additional argument of \code{boot.se}, which is an object returned by \code{\link{boot_lucid}}
#'
#' @return A list, containing the extracted key parameters from the LUCID model that can be used to print the summary table
#' @export
summary.lucid_serial <- function(object, ...){
    args <- list(...)
    auto_print <- if (is.null(args$auto_print)) TRUE else isTRUE(args$auto_print)
    stage_boot <- NULL
    if (!is.null(args$boot.se)) {
      # Preferred format from boot_lucid(lucid_model = "serial")
      if (is.list(args$boot.se) && !is.null(args$boot.se$stage)) {
        stage_boot <- args$boot.se$stage
      } else if (is.list(args$boot.se) && length(args$boot.se) == length(object$submodel)) {
        # Backward-compatible stage-wise list input.
        stage_boot <- args$boot.se
      }
    }

    submodels <- object$submodel
    n_submodels <- length(submodels)
    stage_summary <- vector(mode = "list", length = n_submodels)
    stage_type <- character(n_submodels)
    stage_K <- vector(mode = "list", length = n_submodels)
    transition_labels <- vector(mode = "list", length = n_submodels)
    transition_prev_stage_type <- rep(NA_character_, n_submodels)

    for (i in seq_len(n_submodels)) {
      boot_i <- NULL
      if (!is.null(stage_boot) && length(stage_boot) >= i) {
        boot_i <- stage_boot[[i]]
      }
      stage_summary[[i]] <- summary_lucid(submodels[[i]], boot.se = boot_i)
      stage_type[i] <- if (inherits(submodels[[i]], "early_lucid")) {
        "early"
      } else if (inherits(submodels[[i]], "lucid_parallel")) {
        "parallel"
      } else {
        class(submodels[[i]])[1]
      }
      stage_K[[i]] <- submodels[[i]]$K
    }

    build_transition_labels <- function(prev_model, prev_stage_idx) {
      if (inherits(prev_model, "early_lucid")) {
        k_prev <- as.integer(prev_model$K)
        if (length(k_prev) == 0 || is.na(k_prev) || k_prev <= 1) {
          return(character(0))
        }
        return(paste0("Stage", prev_stage_idx, ".cluster", seq.int(2, k_prev)))
      }
      if (inherits(prev_model, "lucid_parallel")) {
        k_prev <- as.integer(prev_model$K)
        out <- character(0)
        for (layer_idx in seq_along(k_prev)) {
          if (!is.na(k_prev[layer_idx]) && k_prev[layer_idx] > 1) {
            out <- c(out, paste0("Stage", prev_stage_idx,
                                 ".Layer", layer_idx,
                                 ".cluster", seq.int(2, k_prev[layer_idx])))
          }
        }
        return(out)
      }
      character(0)
    }

    transition_labels[[1]] <- character(0)
    for (i in seq.int(2, n_submodels)) {
      transition_labels[[i]] <- build_transition_labels(submodels[[i - 1]], i - 1)
      transition_prev_stage_type[i] <- stage_type[i - 1]
    }

    BIC <- cal_bic_serial(object)
    loglik <- cal_loglik_serial(object)
    reg_summary <- NULL
    if(!is.null(object$Rho)) {
      reg_summary <- list(
        Rho_G = object$Rho$Rho_G,
        Rho_Z_Mu = object$Rho$Rho_Z_Mu,
        Rho_Z_Cov = object$Rho$Rho_Z_Cov
      )
    }

    results <- list(
      model_info = list(
        family = object$family,
        n_observations = object$N,
        n_stages = n_submodels,
        stage_type = stage_type,
        stage_K = stage_K
      ),
      model_fit = list(
        BIC = BIC,
        loglik = loglik
      ),
      regularization = reg_summary,
      missing_data = object$missing_summary,
      stage_summary = stage_summary,
      transition = list(
        labels = transition_labels,
        prev_stage_type = transition_prev_stage_type
      ),
      boot.se = args$boot.se,
      # Keep legacy fields for compatibility with existing downstream code.
      summary.list = stage_summary,
      BIC = BIC,
      loglik = loglik
    )

    class(results) <- "sumlucid_serial"
    if (auto_print) {
      print(results)
    }
    invisible(results)
}


summary_lucid_auxi <- function(object, boot.se = NULL){
  if (inherits(object, "early_lucid")){
    s1 <- object$select$selectG
    s2 <- object$select$selectZ
    nG <- sum(s1)
    nZ <- sum(s2)
    K <- object$K
    gamma <- object$res_Gamma$beta
    #obtain number of parameters
    if(is_normal_family(object$family)){
      nY <- length(object$res_Gamma$beta) + length(object$res_Gamma$sigma)
    }
    if(is_binary_family(object$family)){
      nY <- length(object$res_Gamma$beta)
    }
    npars <- (nG + 1) * (K - 1) + (nZ * K + nZ^2 * K) + nY
    BIC <- -2 * object$likelihood + npars * log(nrow(object$inclusion.p))
    results <- list(beta = object$res_Beta[, c(TRUE, s1)],
                    mu = object$res_Mu[, s2],
                    gamma = object$res_Gamma,
                    family = object$family,
                    K = K,
                    BIC = BIC,
                    loglik = object$likelihood,
                    boot.se = boot.se)
    class(results) <- "sumlucid_early"
    return(results)
  }
  if (inherits(object, "lucid_parallel")){
    ##not having regularity yet, to be added
    s1 <- object$select$selectG
    s2 <- object$select$selectZ
    nG <- sum(s1)
    nZ <- sapply(s2,sum)
    K <- object$K
    #obtain number of parameters
    if(to_parallel_family(object$family) == "gaussian"){
      nY <- length(object$res_Gamma$Gamma$mu) + length(object$res_Gamma$Gamma$sd)
    }
    if(to_parallel_family(object$family) == "binomial"){
      #binary summary res_Gamma$Gamma$mu is unclear, use object$res_Gamma$fit$coefficients instead
      nY <- length(object$res_Gamma$fit$coefficients)
    }
    #initiate number of parameters
    npars = 0
    #compute number of parameters for G to X association
    for (i in 1:length(K)){
      npars_new = (nG + 1) * (K[i] - 1)
      npars = npars + npars_new
    }
    #compute number of parameters for X to Z association and add
    for (i in 1:length(K)){
      npars_new = (nZ[i] * K[i] + nZ[i] * nZ[i] * K[i])
      npars = npars + npars_new
    }
    #compute number of parameters for X to Y association and add
    npars = npars + nY

    BIC <- -2 * object$likelihood + npars * log(nrow(object$inclusion.p[[1]]))

    results <- list(beta = object$res_Beta,
                    mu = object$res_Mu,
                    Gamma = object$res_Gamma,
                    family = object$family,
                    K = K,
                    BIC = BIC,
                    loglik = object$likelihood,
                    #BOOT.SE IS NULL FOR NOW
                    boot.se = NULL
    )
    class(results) <- "sumlucid_parallel"
    return(results)

  }
}

#' @title Print the output of LUCID in a nicer table
#'
#' @param x An object returned by \code{summary}
#' @param ... Other parameters to be passed to \code{print.sumlucid_serial}
#' @return Prints a structured model summary, including model specification, missing-data
#' profile, feature-selection overview, model fit statistics, regularization settings,
#' and detailed parameter estimates. If \code{boot.se} is provided in \code{summary()},
#' bootstrap CI tables are shown for sections (1) Y, (2) Z, and (3) E.
#' @export 
#' @examples
#' \donttest{
#' # use simulated data
#' G <- sim_data$G[1:300, , drop = FALSE]
#' Z <- sim_data$Z[1:300, , drop = FALSE]
#' Y_normal <- sim_data$Y_normal[1:300]
#'
#' # fit lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", family = "normal", K = 2,
#' seed = 1008)
#'
#' # conduct bootstrap resampling
#' boot1 <- suppressWarnings(
#'   boot_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", model = fit1, R = 5)
#' )
#'
#' # print the summary of the lucid model in a table
#' temp <- summary(fit1)
#' print(temp)
#' }
print.sumlucid_early <- function(x, ...) {
  serial_ctx <- x$serial_context
  in_serial <- is.list(serial_ctx) && isTRUE(serial_ctx$serial_stage)
  is_last_stage <- !(in_serial && !isTRUE(serial_ctx$is_last_stage))
  transition_labels <- if (in_serial) serial_ctx$transition_labels else character(0)

  cat("\n====================================================\n")
  cat("LUCID Early Integration: Model Summary\n")
  cat("====================================================\n\n")

  cat("Model specification\n")
  cat(sprintf("  Family                 : %s\n", x$model_info$family))
  cat(sprintf("  Number of observations : %d\n", x$model_info$n_observations))
  cat(sprintf("  Number of clusters (K) : %d\n", x$model_info$K))
  cat("\n")

  if(!is.null(x$missing_data)) {
    md <- x$missing_data
    cat("Missing-data profile\n")
    cat(sprintf("  Listwise missing rows  : %d / %d (%.1f%%)\n",
                md$listwise_rows, md$n_rows, 100 * md$prop_listwise_rows))
    cat(sprintf("  Sporadic missing rows  : %d / %d (%.1f%%)\n",
                md$sporadic_rows, md$n_rows, 100 * md$prop_sporadic_rows))
    cat(sprintf("  Missing cells total    : %d / %d (%.1f%%)\n",
                md$total_missing_cells, md$n_rows * md$n_features,
                100 * md$prop_total_missing_cells))
    cat("\n")
  }

  cat("Feature selection overview\n")
  if(!is.null(x$feature_selection$G)) {
    tab_g <- as.integer(table(x$feature_selection$G$Selected))
    names(tab_g) <- names(table(x$feature_selection$G$Selected))
    g_selected <- ifelse("TRUE" %in% names(tab_g), tab_g["TRUE"], 0)
    g_total <- nrow(x$feature_selection$G)
    g_prop <- ifelse(g_total > 0, g_selected / g_total, 0)
    cat(sprintf("  G features selected    : %d / %d (%.1f%%)\n",
                g_selected, g_total, 100 * g_prop))
  }
  if(!is.null(x$feature_selection$Z)) {
    tab_z <- as.integer(table(x$feature_selection$Z$Selected))
    names(tab_z) <- names(table(x$feature_selection$Z$Selected))
    z_selected <- ifelse("TRUE" %in% names(tab_z), tab_z["TRUE"], 0)
    z_total <- nrow(x$feature_selection$Z)
    z_prop <- ifelse(z_total > 0, z_selected / z_total, 0)
    cat(sprintf("  Z features selected    : %d / %d (%.1f%%)\n",
                z_selected, z_total, 100 * z_prop))
  }
  cat("\n")

  cat("Model fit statistics\n")
  cat(sprintf("  Log-likelihood         : %.2f\n", x$model_fit$loglik))
  cat(sprintf("  BIC                    : %.2f\n", x$model_fit$BIC))
  cat(sprintf("  Number of parameters   : %d\n", x$model_fit$n_parameters))

  if(!is.null(x$regularization)) {
    cat("\nRegularization\n")
    cat(sprintf("  Rho_G                  : %.3f\n", x$regularization$Rho_G))
    cat(sprintf("  Rho_Z_Mu               : %.3f\n", x$regularization$Rho_Z_Mu))
    cat(sprintf("  Rho_Z_Cov              : %.3f\n", x$regularization$Rho_Z_Cov))
  }

  cat("\nDetailed parameter estimates\n")
  K <- x$model_info$K
  beta <- as.data.frame(x$parameters$beta)
  dim1 <- ncol(beta) - 1
  z.mean <- as.data.frame(t(x$parameters$mu))
  
  y <- switch(x$model_info$family, 
              normal = f.normal.early,
              binary = f.binary.early)
  if (is_last_stage) {
    y(x$parameters$gamma, K, se = x$boot.se$gamma)
    cat("\n")
  }
  
  z_idx <- if (is_last_stage) 2 else 1
  cat(sprintf("(%d) Z: mean of omics data for each latent cluster \n", z_idx))
  if(is.null(x$boot.se)){
    colnames(z.mean) <- paste0("mu_cluster", 1:K)
    print(z.mean)
  } else{
    print(x$boot.se$mu)
  }
  cat("\n")
  
  e_idx <- z_idx + 1
  if (in_serial && length(transition_labels) > 0) {
    cat(sprintf("(%d) E: intercept and odds ratio of being assigned to each latent cluster for each cluster from previous serial stage\n", e_idx))
  } else {
    cat(sprintf("(%d) E: intercept and odds ratio of being assigned to each latent cluster for each exposure \n", e_idx))
  }
  if(is.null(x$boot.se)){
    beta_mat <- as.matrix(beta)
    if(is.null(dim(beta_mat))) {
      beta_mat <- matrix(beta_mat, nrow = 1)
    }
    if(ncol(beta_mat) == 0) {
      cat("No E-coefficient estimates are available.\n")
    } else {
      if(nrow(beta_mat) == K) {
        beta_use <- beta_mat[2:K, , drop = FALSE]
        cluster_ids <- 2:K
      } else if(nrow(beta_mat) == (K - 1)) {
        beta_use <- beta_mat
        cluster_ids <- 2:K
      } else {
        beta_use <- beta_mat
        cluster_ids <- seq_len(nrow(beta_use))
      }

      if(nrow(beta_use) == 0) {
        cat("No non-reference cluster coefficient is available.\n")
      } else {
        n_features <- max(0, ncol(beta_use) - 1)
        mapped_names <- character(0)
        if (length(transition_labels) > 0 && n_features > 0) {
          mapped_names <- transition_labels[seq_len(min(length(transition_labels), n_features))]
          if (length(mapped_names) < n_features) {
            mapped_names <- c(mapped_names, paste0("PrevStageCluster", seq.int(length(mapped_names) + 1, n_features)))
          }
        }
        if (length(mapped_names) == n_features && n_features > 0) {
          feature_names <- c("(Intercept)", mapped_names)
        } else {
          feature_names <- colnames(beta_use)
          if(is.null(feature_names)) {
            feature_names <- c("(Intercept)",
                               paste0("G", seq_len(n_features)))
          }
          feature_names[1] <- "(Intercept)"
        }

        g.or <- do.call(rbind, lapply(seq_along(cluster_ids), function(r) {
          data.frame(
            feature = feature_names,
            cluster = paste0("cluster", cluster_ids[r]),
            beta = as.numeric(beta_use[r, ]),
            OR = exp(as.numeric(beta_use[r, ])),
            stringsAsFactors = FALSE
          )
        }))
        rownames(g.or) <- paste0(g.or$feature, ".", g.or$cluster)
        print(g.or[, c("beta", "OR"), drop = FALSE])
      }
    }
  } else{
    beta_ci <- as.data.frame(x$boot.se$beta)
    rn <- rownames(beta_ci)
    if(!is.null(rn)) {
      rn <- sub("^intercept(\\.cluster[0-9]+)$", "(Intercept)\\1",
                rn, ignore.case = TRUE)
      if (length(transition_labels) > 0) {
        for (j in seq_along(transition_labels)) {
          rn <- sub(paste0("^G_?", j, "(\\.cluster[0-9]+)$"),
                    paste0(transition_labels[j], "\\1"), rn)
        }
      }
      rownames(beta_ci) <- rn
    }
    print(beta_ci)
  }

  invisible(x)
}


#' @title Print the output of LUCID in a nicer table
#'
#' @param x An object returned by \code{summary}
#' @param ... Other parameters to be passed to \code{print.sumlucid_parallel}
#' @return Prints a structured parallel-model summary, including per-layer missing-data
#' profile, overall and per-layer feature-selection overview, model fit statistics,
#' regularization settings, and detailed parameter estimates for Y, Z, and E.
#' If \code{boot.se} is provided in \code{summary()}, bootstrap CI tables are
#' shown for sections (1) Y, (2) Z, and (3) E.
#' @export 
print.sumlucid_parallel <- function(x, ...) {
  serial_ctx <- x$serial_context
  in_serial <- is.list(serial_ctx) && isTRUE(serial_ctx$serial_stage)
  is_last_stage <- !(in_serial && !isTRUE(serial_ctx$is_last_stage))
  transition_labels <- if (in_serial) serial_ctx$transition_labels else character(0)

  cat("\n====================================================\n")
  cat("LUCID Parallel: Model Summary\n")
  cat("====================================================\n\n")

  cat("Model specification\n")
  cat(sprintf("  Family                 : %s\n", x$model_info$family))
  cat(sprintf("  Number of observations : %d\n", x$model_info$n_observations))
  cat(sprintf("  Clusters per layer     : %s\n", 
              paste(x$model_info$K, collapse = ", ")))
  cat("\n")

  if(!is.null(x$missing_data) && !is.null(x$missing_data$layer_summary)) {
    cat("Missing-data profile by layer\n")
    md <- x$missing_data$layer_summary
    for (i in seq_len(nrow(md))) {
      cat(sprintf("  Layer %d listwise rows : %d / %d (%.1f%%)\n",
                  md$layer[i], md$listwise_rows[i], md$n_rows[i],
                  100 * md$prop_listwise_rows[i]))
      cat(sprintf("  Layer %d sporadic rows : %d / %d (%.1f%%)\n",
                  md$layer[i], md$sporadic_rows[i], md$n_rows[i],
                  100 * md$prop_sporadic_rows[i]))
      cat(sprintf("  Layer %d missing cells : %d / %d (%.1f%%)\n",
                  md$layer[i], md$total_missing_cells[i],
                  md$n_rows[i] * md$n_features[i],
                  100 * md$prop_total_missing_cells[i]))
    }
    cat("\n")
  }

  cat("Feature selection overview\n")
  if(!is.null(x$feature_selection$G)) {
    n_total <- nrow(x$feature_selection$G)
    n_selected <- sum(x$feature_selection$G$Selected)
    p_selected <- ifelse(n_total > 0, n_selected / n_total, 0)
    cat(sprintf("  G features selected    : %d / %d (%.1f%%)\n",
                n_selected, n_total, 100 * p_selected))
  }
  if(!is.null(x$feature_selection$G_layer)) {
    cat("  G features by layer\n")
    for(i in seq_along(x$feature_selection$G_layer)) {
      n_total <- nrow(x$feature_selection$G_layer[[i]])
      n_selected <- sum(x$feature_selection$G_layer[[i]]$Selected)
      p_selected <- ifelse(n_total > 0, n_selected / n_total, 0)
      cat(sprintf("    Layer %d              : %d / %d (%.1f%%)\n",
                  i, n_selected, n_total, 100 * p_selected))
    }
  }
  if(!is.null(x$feature_selection$Z)) {
    cat("  Z features\n")
    for(i in seq_along(x$feature_selection$Z)) {
      n_total <- nrow(x$feature_selection$Z[[i]])
      n_selected <- sum(x$feature_selection$Z[[i]]$Selected)
      n_multi <- sum(x$feature_selection$Z[[i]]$Selected_in_clusters > 1)
      p_selected <- ifelse(n_total > 0, n_selected / n_total, 0)
      cat(sprintf("    Layer %d selected     : %d / %d (%.1f%%)\n",
                  i, n_selected, n_total, 100 * p_selected))
      cat(sprintf("    Layer %d multi-cluster: %d\n", i, n_multi))
    }
  }
  cat("\n")

  cat("Model fit statistics\n")
  cat(sprintf("  Log-likelihood         : %.2f\n", x$model_fit$loglik))
  cat(sprintf("  BIC                    : %.2f\n", x$model_fit$BIC))
  cat(sprintf("  Number of parameters   : %d\n", x$model_fit$n_parameters))

  if(!is.null(x$regularization)) {
    cat("\nRegularization\n")
    cat(sprintf("  Rho_G                  : %.3f\n", x$regularization$Rho_G))
    cat(sprintf("  Rho_Z_Mu               : %.3f\n", x$regularization$Rho_Z_Mu))
    cat(sprintf("  Rho_Z_Cov              : %.3f\n", x$regularization$Rho_Z_Cov))
  }

  cat("\nDetailed parameter estimates\n")

  # (1) Y model
  se_gamma <- NULL
  if(!is.null(x$boot.se) && !is.null(x$boot.se$gamma)) {
    se_gamma <- x$boot.se$gamma
  }
  y <- switch(to_parallel_family(x$model_info$family),
              gaussian = f.normal.parallel,
              binomial = f.binary.parallel)
  if (is_last_stage) {
    y(x$parameters$Gamma, x$model_info$K, se = se_gamma)
    cat("\n")
  }

  # (2) Z means by layer
  z_idx <- if (is_last_stage) 2 else 1
  cat(sprintf("(%d) Z: mean of omics data for each latent cluster of each layer \n", z_idx))
  se_mu <- NULL
  if(!is.null(x$boot.se) && !is.null(x$boot.se$mu)) {
    se_mu <- x$boot.se$mu
  }
  for(i in seq_along(x$model_info$K)) {
    cat(sprintf("Layer %d\n\n", i))

    if(!is.null(se_mu) && is.list(se_mu) && length(se_mu) >= i && !is.null(se_mu[[i]])) {
      print(se_mu[[i]])
      cat("\n")
      next
    }

    mu_i <- x$parameters$mu[[i]]
    K_i <- x$model_info$K[i]

    if(is.null(dim(mu_i))) {
      mu_i <- matrix(mu_i, nrow = K_i)
    }

    z_mean <- if(nrow(mu_i) == K_i) t(mu_i) else mu_i

    # If selection info is available, display selected Z features only.
    if(!is.null(x$feature_selection$Z) && length(x$feature_selection$Z) >= i) {
      z_sel <- x$feature_selection$Z[[i]]
      if(!is.null(z_sel) && nrow(z_sel) == nrow(z_mean)) {
        rownames(z_mean) <- z_sel$Feature
        keep <- z_sel$Selected
        if(any(keep)) {
          z_mean <- z_mean[keep, , drop = FALSE]
        }
      }
    }

    colnames(z_mean) <- paste0("mu_cluster", seq_len(ncol(z_mean)))
    print(as.data.frame(z_mean))
    cat("\n")
  }

  # (3) G effects by layer
  e_idx <- z_idx + 1
  if (in_serial && length(transition_labels) > 0) {
    cat(sprintf("(%d) E: intercept and odds ratio of being assigned to each latent cluster for each cluster from previous serial stage (for each layer)\n", e_idx))
  } else {
    cat(sprintf("(%d) E: intercept and odds ratio of being assigned to each latent cluster for each exposure for each layer \n", e_idx))
  }
  se_beta <- NULL
  if(!is.null(x$boot.se) && !is.null(x$boot.se$beta)) {
    se_beta <- x$boot.se$beta
  }
  beta_list <- x$parameters$beta$Beta
  if(is.null(beta_list) || length(beta_list) == 0) {
    cat("No G-coefficient estimates are available.\n")
  } else {
    for(i in seq_along(beta_list)) {
      cat(sprintf("Layer %d\n\n", i))

      beta_i <- beta_list[[i]]
      if(is.null(dim(beta_i))) {
        beta_i <- matrix(beta_i, nrow = 1)
      }
      K_i <- x$model_info$K[i]
      if(ncol(beta_i) <= 1) {
        cat("No exposure coefficient available in this layer.\n\n")
        next
      }

      # Keep only non-reference cluster effects where identifiable.
      if(nrow(beta_i) == (K_i - 1)) {
        beta_use <- beta_i
        cluster_ids <- 2:K_i
      } else if(nrow(beta_i) == K_i) {
        beta_use <- beta_i[2:K_i, , drop = FALSE]
        cluster_ids <- 2:K_i
      } else {
        beta_use <- beta_i
        cluster_ids <- seq_len(nrow(beta_use))
      }
      if(nrow(beta_use) == 0) {
        cat("No non-reference cluster coefficient available in this layer.\n\n")
        next
      }

      n_exposure <- if(!is.null(x$feature_selection$G)) nrow(x$feature_selection$G) else (ncol(beta_use) - 1)
      n_exposure <- min(n_exposure, ncol(beta_use) - 1)
      if (length(transition_labels) > 0) {
        g_names <- transition_labels[seq_len(min(length(transition_labels), n_exposure))]
        if (length(g_names) < n_exposure) {
          g_names <- c(g_names, paste0("PrevStageCluster", seq.int(length(g_names) + 1, n_exposure)))
        }
      } else {
        g_names <- colnames(beta_use)[1 + seq_len(n_exposure)]
        if(is.null(g_names)) {
          g_names <- paste0("G", seq_len(n_exposure))
        }
      }

      # If selection info is available, display selected G features only.
      selected_g <- rep(TRUE, length(g_names))
      if(!is.null(x$feature_selection$G_layer) && length(x$feature_selection$G_layer) >= i) {
        g_sel <- x$feature_selection$G_layer[[i]]
        if(!is.null(g_sel) && nrow(g_sel) == length(g_names)) {
          selected_g <- g_sel$Selected
        }
      } else if(!is.null(x$feature_selection$G) && nrow(x$feature_selection$G) == length(g_names)) {
        selected_g <- x$feature_selection$G$Selected
      }

      selected_features <- c("(Intercept)", g_names[selected_g])
      if(length(selected_features) == 1) {
        cat("No exposure is selected in this layer under current Rho_G; showing intercept only.\n\n")
      }

      if(!is.null(se_beta) && is.list(se_beta) && length(se_beta) >= i && !is.null(se_beta[[i]])) {
        beta_ci <- as.data.frame(se_beta[[i]])
        rn <- rownames(beta_ci)
        if(!is.null(rn)) {
          layer_pat <- paste0("^Layer", i, "\\.(.*)\\.cluster([0-9]+)$")
          m <- regexec(layer_pat, rn)
          parsed <- regmatches(rn, m)
          ok <- lengths(parsed) == 3
          if(any(ok)) {
            feature <- vapply(parsed[ok], `[`, character(1), 2)
            cluster <- vapply(parsed[ok], `[`, character(1), 3)
            if (length(transition_labels) > 0) {
              for (j in seq_along(transition_labels)) {
                feature <- sub(paste0("^G_?", j, "$"), transition_labels[j], feature)
              }
            }
            beta_ci <- beta_ci[ok, , drop = FALSE]
            rownames(beta_ci) <- paste0(feature, ".cluster", cluster)
            keep <- feature %in% selected_features
            beta_ci <- beta_ci[keep, , drop = FALSE]
          }
        }
        print(beta_ci)
        cat("\n")
        next
      }

      coef_cols <- c(1, 1 + which(selected_g))
      coef_mat <- beta_use[, coef_cols, drop = FALSE]
      if(is.null(dim(coef_mat))) {
        coef_mat <- matrix(coef_mat, nrow = 1)
      }

      out <- do.call(rbind, lapply(seq_along(cluster_ids), function(r) {
        data.frame(
          feature = selected_features,
          cluster = paste0("cluster", cluster_ids[r]),
          beta = as.numeric(coef_mat[r, ]),
          OR = exp(as.numeric(coef_mat[r, ])),
          stringsAsFactors = FALSE
        )
      }))
      rownames(out) <- paste0(out$feature, ".", out$cluster)
      print(out[, c("beta", "OR"), drop = FALSE])
      cat("\n")
    }
  }

  invisible(x)
}


#' @title Print the output of LUCID in a nicer table
#'
#' @param x An object returned by \code{summary}
#' @param ... Other parameters to be passed to \code{print.sumlucid_serial}
#' @return A nice table/several nice tables of the summary of the LUCID model
#' @export 
print.sumlucid_serial <- function(x, ...){
  if (!inherits(x, "sumlucid_serial")) {
    stop("this function only prints summary of lucid in serial model")
  }

  cat("\n====================================================\n")
  cat("LUCID Serial: Model Summary\n")
  cat("====================================================\n\n")

  n_stages <- if(!is.null(x$model_info$n_stages)) x$model_info$n_stages else length(x$summary.list)
  cat("Model specification\n")
  if(!is.null(x$model_info$family)) {
    cat(sprintf("  Family                 : %s\n", x$model_info$family))
  }
  if(!is.null(x$model_info$n_observations)) {
    cat(sprintf("  Number of observations : %d\n", x$model_info$n_observations))
  }
  cat(sprintf("  Number of stages       : %d\n", n_stages))
  if(!is.null(x$model_info$stage_type) && !is.null(x$model_info$stage_K)) {
    for(i in seq_len(n_stages)) {
      k_i <- x$model_info$stage_K[[i]]
      k_txt <- if(is.list(k_i)) {
        paste(vapply(k_i, function(v) paste(v, collapse = ","), character(1)),
              collapse = " | ")
      } else {
        paste(k_i, collapse = ",")
      }
      cat(sprintf("  Stage %d               : %s (K = %s)\n",
                  i, x$model_info$stage_type[i], k_txt))
    }
  }
  cat("\n")

  if(!is.null(x$missing_data) && !is.null(x$missing_data$stage)) {
    cat("Missing-data profile by stage\n")
    md_stage <- x$missing_data$stage
    for(i in seq_along(md_stage)) {
      cat(sprintf("  Stage %d\n", i))
      md_i <- md_stage[[i]]
      if(is.null(md_i)) {
        cat("    Missing summary unavailable\n")
      } else if(!is.null(md_i$layer_summary)) {
        for(j in seq_len(nrow(md_i$layer_summary))) {
          cat(sprintf("    Layer %d listwise/sporadic rows : %d / %d\n",
                      md_i$layer_summary$layer[j],
                      md_i$layer_summary$listwise_rows[j],
                      md_i$layer_summary$sporadic_rows[j]))
        }
      } else if(!is.null(md_i$listwise_rows) && !is.null(md_i$sporadic_rows)) {
        cat(sprintf("    Listwise rows         : %d\n", md_i$listwise_rows))
        cat(sprintf("    Sporadic rows         : %d\n", md_i$sporadic_rows))
      } else {
        cat("    Missing summary unavailable\n")
      }
    }
    cat("\n")
  }

  if(!is.null(x$model_fit)) {
    cat("Model fit statistics\n")
    cat(sprintf("  Log-likelihood         : %.2f\n", x$model_fit$loglik))
    cat(sprintf("  BIC                    : %.2f\n", x$model_fit$BIC))
  } else {
    cat("Model fit statistics\n")
    cat(sprintf("  Log-likelihood         : %.2f\n", x$loglik))
    cat(sprintf("  BIC                    : %.2f\n", x$BIC))
  }

  if(!is.null(x$regularization)) {
    cat("\nRegularization\n")
    cat(sprintf("  Rho_G                  : %.3f\n", x$regularization$Rho_G))
    cat(sprintf("  Rho_Z_Mu               : %.3f\n", x$regularization$Rho_Z_Mu))
    cat(sprintf("  Rho_Z_Cov              : %.3f\n", x$regularization$Rho_Z_Cov))
  }
  cat("\n")

  cat("Stage-wise detailed parameter estimates\n")
  sum_list <- if(!is.null(x$stage_summary)) x$stage_summary else x$summary.list
  for (i in seq_along(sum_list)) {
    stage_label <- if(!is.null(x$model_info$stage_type)) x$model_info$stage_type[i] else class(sum_list[[i]])[1]
    cat(sprintf("\n--- Stage %d (%s) ---\n", i, stage_label))
    stage_obj <- sum_list[[i]]
    stage_obj$serial_context <- list(
      serial_stage = TRUE,
      is_last_stage = (i == n_stages),
      stage_index = i,
      previous_stage_type = if (i > 1 && !is.null(x$transition$prev_stage_type)) x$transition$prev_stage_type[i] else NULL,
      transition_labels = if (i > 1 && !is.null(x$transition$labels)) x$transition$labels[[i]] else character(0)
    )
    print(stage_obj)
  }

  invisible(x)
}



#' @export
print.sumlucid_auxi <- function(x, ...){
  if(inherits(x, "sumlucid_early")){
  K <- x$K
  beta <- as.data.frame(x$beta)
  dim1 <- ncol(beta) - 1
  z.mean <- as.data.frame(t(x$mu))
  cat("----------Summary of the LUCID Early Integration model---------- \n \n")
  cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")
  y <- switch(x$family, normal = f.normal.early,
              binary = f.binary.early)
  y(x$gamma, K, se = x$boot.se$gamma)
  cat("\n")
  cat("(2) Z: mean of omics data for each latent cluster \n")
  if(is.null(x$boot.se)){
    colnames(z.mean) <- paste0("mu_cluster", 1:K)
    print(z.mean)
  } else{
    print(x$boot.se$mu)
  }
  cat("\n")
  cat("(3) E: intercept and odds ratio of being assigned to each latent cluster for each exposure \n")
  if(is.null(x$boot.se)){
    beta_mat <- as.matrix(beta)
    if(is.null(dim(beta_mat))) {
      beta_mat <- matrix(beta_mat, nrow = 1)
    }
    if(ncol(beta_mat) == 0) {
      cat("No E-coefficient estimates are available.\n")
    } else {
      if(nrow(beta_mat) == K) {
        beta_use <- beta_mat[2:K, , drop = FALSE]
        cluster_ids <- 2:K
      } else if(nrow(beta_mat) == (K - 1)) {
        beta_use <- beta_mat
        cluster_ids <- 2:K
      } else {
        beta_use <- beta_mat
        cluster_ids <- seq_len(nrow(beta_use))
      }
      if(nrow(beta_use) == 0) {
        cat("No non-reference cluster coefficient is available.\n")
      } else {
        feature_names <- colnames(beta_use)
        if(is.null(feature_names)) {
          feature_names <- c("(Intercept)",
                             paste0("G", seq_len(max(0, ncol(beta_use) - 1))))
        }
        feature_names[1] <- "(Intercept)"
        g.or <- do.call(rbind, lapply(seq_along(cluster_ids), function(r) {
          data.frame(
            feature = feature_names,
            cluster = paste0("cluster", cluster_ids[r]),
            beta = as.numeric(beta_use[r, ]),
            OR = exp(as.numeric(beta_use[r, ])),
            stringsAsFactors = FALSE
          )
        }))
        rownames(g.or) <- paste0(g.or$feature, ".", g.or$cluster)
        print(g.or[, c("beta", "OR"), drop = FALSE])
      }
    }
  } else{
    beta_ci <- as.data.frame(x$boot.se$beta)
    rn <- rownames(beta_ci)
    if(!is.null(rn)) {
      rn <- sub("^intercept(\\.cluster[0-9]+)$", "(Intercept)\\1",
                rn, ignore.case = TRUE)
      rownames(beta_ci) <- rn
    }
    print(beta_ci)
      }
  }
  if(inherits(x, "sumlucid_parallel")){
    K <- x$K
    #list of betas for each layer
    beta <- x$beta$Beta
    dim1 <- ncol(x$beta$Beta[[1]]) - 1
    #list of Z means for each layer, transposed
    z.mean <- x$mu
    cat("----------Summary of the LUCID in Parallel model---------- \n \n")
    cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")
    y <- switch(to_parallel_family(x$family), gaussian = f.normal.parallel,
                binomial = f.binary.parallel)
    y(x$Gamma, K, se = x$boot.se$gamma)
    cat("\n")
    cat("(2) Z: mean of omics data for each latent cluster of each layer \n")
    if(is.null(x$boot.se)){
      for (i in 1:length(K)){
        cat("Layer ",i, "\n \n")
        colnames(z.mean[[i]]) <- paste0("mu_cluster", 1:K[i])
        print(z.mean[[i]])
        cat("\n \n")
      }

    }else{
      print(x$boot.se$mu)
    }
    cat("\n")
    cat("(3) E: odds ratio of being assigned to each latent cluster for each exposure for each layer \n")
    if(is.null(beta)) {
      cat("no exposure is selected given current penalty Rho_G, please use a smaller penalty")
    } else {
      for (i in 1:length(K)){
      cat("Layer ",i, "\n \n")
      dd <- as.matrix(as.data.frame(beta[[i]][, -1]))
      g.or <- data.frame(beta = unlist(split(dd, row(dd))))

      rownames(g.or) <- paste0(colnames(beta[[i]])[-1], ".cluster", sapply(2:K[i], function(x) return(rep(x, dim1))))


      if(is.null(x$boot.se)){
        g.or$OR <- exp(g.or$beta)
        print(g.or)
      } else{
        print(x$boot.se$beta)
      }
      cat("\n \n")
      }
    }
  }
}


# summarize output of normal outcome for early
f.normal.early <- function(x, K, se){

  cat("(1) Y (continuous outcome): effect size of Y for each latent cluster (and effect of covariates if included) \n")

  normalize_early_gamma_names <- function(nm, K, n_par) {
    if(is.null(nm)) {
      if(n_par >= K) {
        nm <- c(paste0("cluster", seq_len(K)),
                if(n_par > K) paste0("CoY", seq_len(n_par - K)))
      } else {
        nm <- paste0("param", seq_len(n_par))
      }
    }
    nm <- gsub("^LC([0-9]+)$", "cluster\\1", nm)
    nm[nm %in% c("Intercept")] <- "(Intercept)"
    nm
  }

  if(!is.null(se)){
    y <- as.data.frame(se)
    rownames(y) <- normalize_early_gamma_names(rownames(y), K, nrow(y))
    if(ncol(y) >= 1 && identical(colnames(y)[1], "estimate")) {
      colnames(y)[1] <- "Gamma"
    }
  } else {
    beta <- as.numeric(x$beta)
    rn <- names(x$beta)
    if(is.null(rn) && !is.null(dim(x$beta))) {
      rn <- rownames(x$beta)
    }
    rn <- normalize_early_gamma_names(rn, K, length(beta))
    y <- data.frame(Gamma = beta, row.names = rn, check.names = FALSE)
  }
  print(y)
}


# summarize output of binary outcome for early
f.binary.early <- function(x, K, se){
  cat("(1) Y (binary outcome): log odds of Y for cluster 1 (reference) and log OR for rest cluster (and log OR of covariate if included)\n")
  gamma <- as.data.frame(x$beta)
  gamma[1,] = 0
  colnames(gamma) <- "gamma"
  if(is.null(se)){
    gamma$`exp(gamma)` <- exp(gamma$gamma)
  } else{
    gamma <- cbind(gamma, se[, -1])
  }
  print(gamma)
}

normalize_parallel_y_names <- function(nm, K){
  if(is.null(nm)) {
    return(paste0("param", seq_len(sum(pmax(K - 1, 0)) + 1)))
  }

  nm <- gsub("^Y\\.", "", nm)
  total_nonref <- sum(pmax(K - 1, 0))

  map_nonref_index <- function(idx) {
    idx <- as.integer(idx)
    rem <- idx
    for(i in seq_along(K)) {
      n_i <- K[i] - 1
      if(rem <= n_i) {
        return(paste0("cluster", rem + 1, "Layer", i))
      }
      rem <- rem - n_i
    }
    paste0("LC", idx)
  }

  vapply(nm, function(s) {
    if(s %in% c("(Intercept)", "Intercept")) {
      return("(Intercept)")
    }
    if(grepl("^(LC|r_fit)[0-9]+$", s)) {
      idx <- as.integer(sub("^(?:LC|r_fit)([0-9]+)$", "\\1", s))
      return(map_nonref_index(idx))
    }
    if(grepl("^gamma[0-9]+$", s)) {
      idx <- as.integer(sub("^gamma([0-9]+)$", "\\1", s))
      if(idx == 1L) {
        return("(Intercept)")
      }
      if((idx - 1L) <= total_nonref) {
        return(map_nonref_index(idx - 1L))
      }
      return(paste0("cov", idx - total_nonref - 1L))
    }
    s
  }, character(1))
}

# summarize output of normal outcome for parallel
f.normal.parallel <- function(x, K, se){

  cat("(1) Y (continuous outcome): intercept, effects of each non-reference latent cluster for each layer of Y (and effect of covariates if included) \n")

  if(!is.null(se)){
    gamma <- as.data.frame(se)
    if(ncol(gamma) >= 1 && identical(colnames(gamma)[1], "estimate")) {
      colnames(gamma)[1] <- "Gamma"
    }
    rownames(gamma) <- normalize_parallel_y_names(rownames(gamma), K)
  } else {
    coef_vec <- coef(x$fit)
    gamma <- data.frame(Gamma = as.numeric(coef_vec), check.names = FALSE)
    rownames(gamma) <- normalize_parallel_y_names(names(coef_vec), K)
  }
  print(gamma)
}


# summarize output of binary outcome for parallel
f.binary.parallel <- function(x, K, se){
  cat("(1) Y (binary outcome): intercept and log OR for non-reference clusters for each layer (and log OR of covariate if included)\n")
  coef_vec <- coef(x$fit)
  gamma <- data.frame(gamma = as.numeric(coef_vec), check.names = FALSE)
  rownames(gamma) <- normalize_parallel_y_names(names(coef_vec), K)
  if(is.null(se)){
    gamma$`exp(gamma)` <- exp(gamma$gamma)
  } else{
    se_df <- as.data.frame(se)
    rownames(se_df) <- normalize_parallel_y_names(rownames(se_df), K)
    idx <- match(rownames(gamma), rownames(se_df))
    se_df <- se_df[idx, , drop = FALSE]
    gamma <- cbind(gamma, se_df[, -1, drop = FALSE])
  }
  print(gamma)
}

##########functions for LUCID in parallel##########




# rearrange cluster order
#
# for continuous outcome - use the cluster combination corresponding to smallest
# mean as the reference cluster
get_ref_cluster <- function(Delta) {
  K <- Delta$K
  mu <- Delta$mu
  mu_matrix <- vec_to_array(K = K, mu = mu)
  ref_index <- which(mu_matrix == min(mu_matrix))
  ref <- arrayInd(ref_index, .dim = K)
  return(ref)
}


# re-arrange parameters for Delta
reorder_Delta <- function(ref, Delta) {
  K <- Delta$K
  mu_matrix <- vec_to_array(K = K, mu = Delta$mu)
  mu <- mu_matrix[ref]

  # if 1 omics layers
  if(length(K) == 1) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1] - mu_matrix[ref[1], 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])


    K_order <- list(K1 = k1_order)
  }

  # if 2 omics layers
  if(length(K) == 2) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1] - mu_matrix[ref[1], 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i] - mu_matrix[1, ref[2]]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order)
  }


  # if 3 omics layers
  if(length(K) == 3) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1] - mu_matrix[ref[1], 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1] - mu_matrix[1, ref[2], 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i] - mu_matrix[1, 1, ref[3]]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])


    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order)
  }


  # if 4 omics layers
  if(length(K) == 4) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1,1] - mu_matrix[ref[1], 1, 1,1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1,1] - mu_matrix[1, ref[2], 1,1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i,1] - mu_matrix[1, 1, ref[3],1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])

    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i] - mu_matrix[1, 1, 1, ref[4]]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order)
  }

  # if 5 omics layers
  if(length(K) == 5) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1, 1, 1] - mu_matrix[ref[1], 1, 1, 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1, 1, 1] - mu_matrix[1, ref[2], 1, 1, 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i, 1, 1] - mu_matrix[1, 1, ref[3], 1, 1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])

    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i, 1] - mu_matrix[1, 1, 1, ref[4], 1]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])

    # reorder K5
    mu_k5 <- rep(0, K[5])
    for(i in 1:K[5]) {
      mu_k5[i] <- mu_matrix[1, 1, 1, 1, i] - mu_matrix[1, 1, 1, 1, ref[5]]
    }
    mu_k5_sort <- sort(mu_k5)
    mu <- c(mu, mu_k5_sort[mu_k5_sort != 0])
    # order of re-arranged cluster for omics 5
    k5 <- order(mu_k5)
    k5_order <- c(ref[5], k5[k5 != ref[5]])


    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order,
                    K5 = k5_order)
  }


  # if 6 omics layers
  if(length(K) == 6) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1, 1, 1, 1] - mu_matrix[ref[1], 1, 1, 1, 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1, 1, 1, 1] - mu_matrix[1, ref[2], 1, 1, 1, 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i, 1, 1, 1] - mu_matrix[1, 1, ref[3], 1, 1, 1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])

    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i, 1, 1] - mu_matrix[1, 1, 1, ref[4], 1, 1]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])

    # reorder K5
    mu_k5 <- rep(0, K[5])
    for(i in 1:K[5]) {
      mu_k5[i] <- mu_matrix[1, 1, 1, 1, i, 1] - mu_matrix[1, 1, 1, 1, ref[5], 1]
    }
    mu_k5_sort <- sort(mu_k5)
    mu <- c(mu, mu_k5_sort[mu_k5_sort != 0])
    # order of re-arranged cluster for omics 5
    k5 <- order(mu_k5)
    k5_order <- c(ref[5], k5[k5 != ref[5]])

    # reorder K6
    mu_k6 <- rep(0, K[6])
    for(i in 1:K[6]) {
      mu_k6[i] <- mu_matrix[1, 1, 1, 1, 1, i] - mu_matrix[1, 1, 1, 1, 1, ref[6]]
    }
    mu_k6_sort <- sort(mu_k6)
    mu <- c(mu, mu_k6_sort[mu_k6_sort != 0])
    # order of re-arranged cluster for omics 6
    k6 <- order(mu_k6)
    k6_order <- c(ref[6], k6[k6 != ref[6]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order,
                    K5 = k5_order,
                    K6 = k6_order)
  }


  Delta$mu <- mu
  return(list(Delta = Delta,
              K_order = K_order))
}



reorder_Mu_Sigma <- function(Mu_Sigma, K_order) {
  for(i in 1:length(K_order)) {
    temp_Mu <- Mu_Sigma$Mu[[i]]
    temp_Sigma <- Mu_Sigma$Sigma[[i]]
    # reorder Mu
    Mu_Sigma$Mu[[i]] <- temp_Mu[, K_order[[i]]]
    Mu_Sigma$Sigma[[i]] <- temp_Sigma[, , K_order[[i]]]
  }
  return(Mu_Sigma)
}



reorder_Beta <- function(Beta, K_order) {
  for(i in 1:length(K_order)) {
    temp_Beta <- Beta[[i]]
    temp_Beta <- rbind(rep(0, ncol(temp_Beta)),
                       temp_Beta)
    temp_Beta_reorder <- temp_Beta[K_order[[i]], ]
    ref <- temp_Beta_reorder[1, ]
    for(j in 1:nrow(temp_Beta_reorder)) {
      temp_Beta_reorder[j, ] <- temp_Beta_reorder[j, ] - ref
    }
    Beta[[i]] <- temp_Beta_reorder[-1, ]
  }

  return(Beta)
}



reorder_z <- function(z, K_order) {
  if(length(K_order) == 2) {
    z <- z[K_order[[1]], K_order[[2]], ]
  }
  return(z)
}


###function to reorder all model parameters###
reorder_lucid <- function(model) {
  ref <- get_ref_cluster(Delta = model$res_Delta$Delta)
  r_Delta <- reorder_Delta(ref = ref,
                           Delta = model$res_Delta$Delta)
  r_Mu_Sigma <- reorder_Mu_Sigma(model$res_Mu_Sigma,
                                 K_order = r_Delta$K_order)
  r_Beta <- reorder_Beta(Beta = model$res_Beta$Beta,
                         K_order = r_Delta$K_order)
  model$res_Delta$Delta <- r_Delta
  model$res_Mu_Sigma$Mu <- r_Mu_Sigma$Mu
  model$res_Mu_Sigma$Sigma <- r_Mu_Sigma$Sigma
  model$res_Beta$Beta <- r_Beta
  model$z <- reorder_z(model$z, K_order = r_Delta$K_order)
  return(model)
}


# function to calculate BIC for LUCID in parallel
cal_bic_parallel <- function(object) {
  ##not having regularity yet, to be added
  s1 <- object$select$selectG
  s2 <- object$select$selectZ
  nG <- if (is.list(s1)) sum(sapply(s1, sum)) else sum(s1)
  nZ <- sapply(s2, function(x) {
    if (is.null(dim(x))) {
      sum(x)
    } else {
      sum(colSums(x) > 0)
    }
  })
  K <- object$K
  #obtain number of parameters
  if(to_parallel_family(object$family) == "gaussian"){
    nY <- length(object$res_Gamma$Gamma$mu) + length(object$res_Gamma$Gamma$sd)
  }
  if(to_parallel_family(object$family) == "binomial"){
    #binary summary res_Gamma$Gamma$mu is unclear, use object$res_Gamma$fit$coefficients instead
    nY <- length(object$res_Gamma$fit$coefficients)
  }
  #initiate number of parameters
  npars = 0
  #compute number of parameters for G to X association
  for (i in 1:length(K)){
    npars_new = (nG + 1) * (K[i] - 1)
    npars = npars + npars_new
  }
  #compute number of parameters for X to Z association and add
  for (i in 1:length(K)){
    npars_new = (nZ[i] * K[i] + nZ[i] * nZ[i] * K[i])
    npars = npars + npars_new
  }
  #compute number of parameters for X to Y association and add
  npars = npars + nY

  BIC <- -2 * object$likelihood + npars * log(nrow(object$inclusion.p[[1]]))

  return(BIC)
}

# function to calculate BIC for LUCID in serial
cal_bic_serial <- function(object) {

  sum_all_sub = summary_lucid_simple(object)
  BIC = 0
  for (i in 1:length(sum_all_sub)){
    BIC_temp = sum_all_sub[[i]]$BIC
    BIC = BIC + BIC_temp
    }
  return(BIC)
}

cal_loglik_serial <- function(object) {

  sum_all_sub = summary_lucid_simple(object)
  loglik = 0
  for (i in 1:length(sum_all_sub)){
    loglik_temp = sum_all_sub[[i]]$loglik
    loglik = loglik + loglik_temp
  }
  return(loglik)
}

#' Summarize results of the early LUCID model
#'
#' @description
#' Assembles the reported quantities for a fitted LUCID model and, by default,
#' prints them. The same components are returned invisibly as a list of class
#' \code{sumlucid_early}, so they can be extracted programmatically rather than
#' parsed from the printed output.
#'
#' Two conventions are worth noting when reading the output. Outcome effects are
#' printed as an intercept -- cluster 1's level -- followed by explicit
#' contrasts of each remaining cluster against it, so the second row is a
#' between-cluster difference and not that cluster's own mean. And the parameter
#' tables are restricted to the features the model retained, so their dimensions
#' match the fit rather than the original input.
#'
#' @param object A LUCID model fitted by \code{\link{estimate_lucid}} or
#'   \code{\link{lucid}}.
#' @param ... Additional arguments. \code{boot.se} accepts an object returned by
#'   \code{\link{boot_lucid}}, whose confidence limits are then shown alongside
#'   the point estimates. \code{auto_print = FALSE} suppresses printing and
#'   returns the summary list only.
#' @return A list of class \code{sumlucid_early} with components:
#'   \describe{
#'     \item{BIC, loglik}{The Bayesian information criterion and the
#'       observed-data log-likelihood at the estimates. Also repeated inside
#'       \code{model_fit}, alongside \code{n_parameters}, the effective
#'       parameter count the BIC charges (Eq 13, reduced per Eq 18 when a
#'       penalty deselected variables).}
#'     \item{model_info}{The outcome \code{family}, the number of clusters
#'       \code{K}, \code{n_observations}, and \code{n_features}, which counts
#'       \emph{retained} exposures and omics features.}
#'     \item{feature_selection}{Which exposures and omics features survived, and
#'       how many were dropped.}
#'     \item{regularization}{The penalties in force, \code{Rho_G},
#'       \code{Rho_Z_Mu} and \code{Rho_Z_Cov}.}
#'     \item{parameters}{Estimates restricted to the retained features:
#'       \code{beta} (exposure-to-cluster, intercept column always kept along
#'       with any covariate columns), \code{mu} (cluster-specific omics means)
#'       and \code{gamma} (cluster-to-outcome, in both absolute and
#'       reference-coded form).}
#'     \item{missing_data}{The fit's \code{missing_summary}; see
#'       \code{\link{estimate_lucid}}.}
#'     \item{boot.se}{The \code{boot.se} argument as supplied, or \code{NULL}.}
#'   }
#'   When \code{boot.se} is supplied, every printed bootstrap CI table
#'   (G-to-X, cluster-to-Y, and cluster-specific omics means) gains a
#'   \code{sig} column: \code{"*"} where the normal-theory confidence
#'   interval excludes 0, \code{""} otherwise.
#' @seealso \code{\link{boot_lucid}} for the confidence limits, and
#'   \code{\link{predict_lucid}} for cluster and outcome prediction.
#' @export
#' @examples 
#' \donttest{
#' # use simulated data (a small subset keeps the example quick)
#' G <- sim_data$G[1:150, , drop = FALSE]
#' Z <- sim_data$Z[1:150, , drop = FALSE]
#' Y_normal <- sim_data$Y_normal[1:150]
#'
#' # fit lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", family = "normal", K = 2,
#' seed = 1008, max_itr = 20, max_tot.itr = 50)
#'
#' # conduct bootstrap resampling
#' boot1 <- suppressWarnings(
#'   boot_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", model = fit1, R = 3)
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

#' Mark rows of a bootstrap CI table whose interval excludes the null value
#'
#' Every bootstrap CI table \code{summary()} prints (G-to-X, cluster-to-Y,
#' and cluster-specific omics means, for both early and parallel, normal and
#' binary outcome) is on the linear-predictor or raw-mean scale, never
#' already exponentiated to an odds ratio, so a single null value, 0, and a
#' single pair of CI columns (\code{gen_ci()}'s
#' \code{norm_lower}/\code{norm_upper}) covers every table the package
#' prints. Adds a \code{sig} column, \code{"*"} where the interval excludes
#' 0, \code{""} otherwise. Also drops the percentile-interval columns
#' (\code{perc_lower}/\code{perc_upper}) from the printed table: \code{sig}
#' is always determined from the normal-theory interval alone, so showing a
#' second, unused pair of bounds next to it only invites the reader to
#' (wrongly) think significance was assessed against them too. The
#' percentile interval is still computed and available on the raw
#' \code{boot_lucid()} object -- only the printed summary table is trimmed.
#'
#' @param df A data frame with \code{lower}/\code{upper} columns from
#'   \code{gen_ci()} (or a subset/reordering of one).
#' @param lower,upper Column names holding the CI bounds to test.
#' @return \code{df} with the percentile-interval columns removed (if
#'   present) and one appended column, \code{sig}.
#' @noRd
.append_significance <- function(df, lower = "norm_lower", upper = "norm_upper") {
  if (is.null(df) || !all(c(lower, upper) %in% colnames(df))) return(df)
  # A CI table is sometimes a plain matrix (e.g. the omics-mean table) rather
  # than a data frame; `[[` on a matrix does not do name-based column lookup
  # ("subscript out of bounds"), so normalize to a data frame first.
  df <- as.data.frame(df)
  df <- df[, setdiff(colnames(df), c("perc_lower", "perc_upper")), drop = FALSE]
  excludes_null <- df[[lower]] > 0 | df[[upper]] < 0
  df$sig <- ifelse(is.na(excludes_null), "", ifelse(excludes_null, "*", ""))
  df
}

#' Number of estimated outcome-model parameters
#'
#' @param object A fitted \code{early_lucid} or \code{lucid_parallel} object.
#' @return \code{0} if \code{object} is unsupervised, otherwise the count of
#'   outcome coefficients (plus one residual-variance parameter for a normal
#'   outcome).
#' @noRd
outcome_npars <- function(object) {
  if (!isTRUE(object$useY)) return(0L)
  if (inherits(object, "early_lucid")) {
    return(length(object$res_Gamma$beta) +
             if (is_normal_family(object$family)) length(object$res_Gamma$sigma) else 0L)
  }
  length(parallel_delta_coef(object$res_Gamma$Gamma)) +
    as.integer(to_parallel_family(object$family) == "gaussian")
}

#' Number of free entries across K symmetric covariance matrices
#'
#' @param n_features Number of features (matrix dimension).
#' @param K Number of clusters.
#' @return \code{K * n_features * (n_features + 1) / 2}.
#' @noRd
covariance_npars <- function(n_features, K) {
  K * n_features * (n_features + 1L) / 2L
}

#' Effective number of parameters for the early model (Eqs 13, 18)
#'
#' Unpenalized: \eqn{D = (P + nCoG + 1)(K - 1) + KM + KM(M + 1)/2 + nY} (Eq
#' 13). Penalized: BICp uses \eqn{D - D_G - D_Z} (Eq 18), where \eqn{D_G}
#' and \eqn{D_Z} count the exposures/omics features whose effects are zero
#' in every cluster -- Eq 18 subtracts one per deselected variable, not that
#' variable's whole mean/covariance block.
#'
#' Eq 13 as printed omits the multinomial-logit intercepts. They are
#' genuinely estimated parameters, so they are counted here, making \code{D}
#' larger than the paper's expression by \eqn{K - 1}.
#'
#' @param object A fitted \code{early_lucid} object.
#' @return The effective parameter count, an integer.
#' @noRd
early_npars <- function(object) {
  selected_g <- object$select$selectG
  selected_z <- object$select$selectZ
  K <- object$K
  n_exposure <- length(selected_g)
  n_cog <- max(0L, length(object$var.names$Gnames) - n_exposure)
  n_z_total <- ncol(object$Z)

  D <- (n_exposure + n_cog + 1L) * (K - 1L) +
    K * n_z_total + covariance_npars(n_z_total, K) + outcome_npars(object)

  DG <- if (is.null(selected_g)) 0L else sum(!selected_g)
  DZ <- if (is.null(selected_z)) 0L else sum(!selected_z)
  D - DG - DZ
}

#' Effective number of parameters for the parallel model (Eqs 13, 18, per layer)
#'
#' Same formula as \code{\link{early_npars}}, applied per layer and summed.
#'
#' @param object A fitted \code{lucid_parallel} object.
#' @return The effective parameter count, an integer.
#' @noRd
parallel_npars <- function(object) {
  K <- object$K
  selected_g <- object$select$selectG_layer
  if (is.null(selected_g)) selected_g <- object$select$selectG
  if (!is.list(selected_g)) selected_g <- rep(list(selected_g), length(K))
  n_exposure <- if (length(selected_g)) length(selected_g[[1]]) else 0L
  n_cog <- max(0L, length(object$var.names$Gnames) - n_exposure)
  # Eq 18 applied per layer: full dimension, minus one per deselected variable.
  DG <- vapply(selected_g, function(x) if (is.null(x)) 0 else sum(!x), numeric(1))
  selected_z <- object$select$selectZ
  n_z_total <- vapply(seq_along(K), function(i) ncol(object$Z[[i]]), numeric(1))
  DZ <- vapply(seq_along(K), function(i) {
    x <- selected_z[[i]]
    if (is.null(x)) return(0)
    sum(!.selectZ_feature_mask(x))
  }, numeric(1))
  D <- sum((n_exposure + n_cog + 1L) * (K - 1L)) +
    sum(K * n_z_total + covariance_npars(n_z_total, K)) + outcome_npars(object)
  D - sum(DG) - sum(DZ)
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
  
  npars <- early_npars(object)
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
      beta = object$res_Beta[, c(TRUE, selectG,
        rep(TRUE, max(0L, ncol(object$res_Beta) - length(selectG) - 1L))), drop = FALSE],
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
#' @return A list of class \code{sumlucid_parallel}, with the same components as
#'   the early-model summary (see \code{\link{summary_lucid}}) but resolved per
#'   omics layer: \code{feature_selection} reports both the overall retained set
#'   and the per-layer sets, \code{parameters} holds \code{mu} by layer and
#'   \code{beta} by layer alongside the single \code{gamma}, and
#'   \code{missing_data} is a list by layer. \code{BIC} and \code{loglik} remain
#'   scalars for the joint model. Returned invisibly when printed.
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
    sapply(selectZ, function(x) sum(.selectZ_feature_mask(x)))
  }
  
  npars <- parallel_npars(object)
  
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
      selected_counts <- if (is.null(dim(selected_i))) {
        as.integer(selected_i)
      } else {
        colSums(selected_i)
      }
      selected_any <- .selectZ_feature_mask(selected_i)
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
    # Top-level BIC/loglik mirror the early-model summary so that
    # summary(fit)$BIC works for every model type.  They previously existed only
    # for "early", so the same expression silently returned NULL for parallel.
    BIC = BIC,
    loglik = object$likelihood,
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
#' @param ... Additional arguments. \code{boot.se} accepts an object returned by
#'   \code{\link{boot_lucid}}, whose confidence limits are then shown alongside
#'   the point estimates. \code{auto_print = FALSE} suppresses printing.
#'
#' @return A list of class \code{sumlucid_serial} with components:
#'   \describe{
#'     \item{BIC, loglik}{Assembled over the whole serial model, and repeated
#'       inside \code{model_fit}.}
#'     \item{model_info}{The outcome \code{family}, \code{n_observations},
#'       \code{n_stages}, \code{stage_type} (whether each stage is an "early" or
#'       a "parallel" sub-model) and \code{stage_K}, the clusters per stage.}
#'     \item{regularization}{The penalties in force, or \code{NULL} if none.}
#'     \item{missing_data}{\code{n_stages} and a per-stage breakdown.}
#'     \item{stage_summary}{The heart of the object: one summary per stage, each
#'       with the same shape as the corresponding early or parallel summary.
#'       Outcome parameters appear on the final stage, the only one the outcome
#'       enters.}
#'     \item{transition}{How consecutive stages connect: \code{labels} names the
#'       previous stage's clusters as they enter the next stage's design, and
#'       \code{prev_stage_type} records that stage's type. The first element is
#'       empty, the first stage having no predecessor.}
#'     \item{boot.se}{The \code{boot.se} argument as supplied, or \code{NULL}.}
#'     \item{summary.list}{A legacy alias of \code{stage_summary}, retained for
#'       downstream code written before the rename. Prefer
#'       \code{stage_summary}.}
#'   }
#'   Returned invisibly when printed.
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
      # see the parallel branch: keep summary(fit)$BIC available for all types
      BIC = BIC,
      loglik = loglik,
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
      # summary.list is a legacy alias of stage_summary, kept for downstream
      # code that predates the rename.  BIC and loglik were repeated here as
      # well, which left the returned list with two entries under each of those
      # names -- `$BIC` resolved to the first, so the copies were unreachable
      # and only made str() output look malformed.  They are already at the top
      # of the list.
      summary.list = stage_summary
    )

    class(results) <- "sumlucid_serial"
    if (auto_print) {
      print(results)
    }
    invisible(results)
}


#' Simplified per-model summary, without printing
#'
#' Routes to \code{\link{summary_lucid_auxi}} directly for early/parallel,
#' or once per submodel for serial. Used by \code{\link{cal_bic_serial}} and
#' \code{\link{cal_loglik_serial}} to aggregate across stages.
#'
#' @param object A fitted LUCID model.
#' @param boot.se Optional \code{boot_lucid()} result to attach.
#' @return For early/parallel, the result of \code{summary_lucid_auxi()}.
#'   For serial, a list of class \code{sumlucid_serial}, one such result per
#'   submodel.
#' @noRd
summary_lucid_simple <- function(object, boot.se = NULL){
  if (inherits(object, "early_lucid") | inherits(object, "lucid_parallel")){
    summary_lucid_auxi(object = object, boot.se = boot.se)
  }
  else if (inherits(object, "lucid_serial")){
    K = object$K
    submodels = object$submodel
    n_submodels = length(submodels)
    summary.list <- vector(mode = "list", length = n_submodels)
    for (i in 1:n_submodels){
      summary.list[[i]] = summary_lucid_auxi(submodels[[i]])
    }

    #return a list of sumlucid object for each submodel
    class(summary.list) <- "sumlucid_serial"
    return(summary.list)
  }
}


#' Assemble the reported quantities for one early/parallel model (or submodel)
#'
#' The shared workhorse behind \code{summary.early_lucid()},
#' \code{summary.lucid_parallel()}, and (per submodel, via
#' \code{\link{summary_lucid_simple}}) \code{summary.lucid_serial()}.
#'
#' @param object A fitted \code{early_lucid} or \code{lucid_parallel} object
#'   (or serial submodel of one of those classes).
#' @param boot.se Optional \code{boot_lucid()} result to attach.
#' @return A list of class \code{sumlucid_early}/\code{sumlucid_parallel};
#'   see \code{\link{summary.early_lucid}}/\code{\link{summary.lucid_parallel}}
#'   for the component-by-component description.
#' @noRd
summary_lucid_auxi <- function(object, boot.se = NULL){
  if (inherits(object, "early_lucid")){
    s1 <- object$select$selectG
    s2 <- object$select$selectZ
    nG <- sum(s1)
    nZ <- sum(s2)
    K <- object$K
    gamma <- object$res_Gamma$beta
    npars <- early_npars(object)
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
    npars <- parallel_npars(object)

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
#' # use simulated data (a small subset keeps the example quick)
#' G <- sim_data$G[1:150, , drop = FALSE]
#' Z <- sim_data$Z[1:150, , drop = FALSE]
#' Y_normal <- sim_data$Y_normal[1:150]
#'
#' # fit lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", family = "normal", K = 2,
#' seed = 1008, max_itr = 20, max_tot.itr = 50)
#'
#' # conduct bootstrap resampling
#' boot1 <- suppressWarnings(
#'   boot_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", model = fit1, R = 2)
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
    print(.append_significance(x$boot.se$mu))
  }
  cat("\n")
  
  e_idx <- z_idx + 1
  .e_source <- if (in_serial && length(transition_labels) > 0) {
    "each cluster from the previous serial stage"
  } else {
    "each exposure"
  }
  if(is.null(x$boot.se)){
    cat(sprintf("(%d) E: intercept and odds ratio of being assigned to each latent cluster for %s\n",
                e_idx, .e_source))
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
    cat(sprintf("(%d) E: cluster-assignment coefficients (log-odds scale) for %s, with bootstrap interval and the odds ratio (exp of each column). `sig` marks intervals excluding the null.\n",
                e_idx, .e_source))
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
    print(.with_or_columns(.append_significance(beta_ci)))
  }

  invisible(x)
}


#' @title Print the output of LUCID in a nicer table
#'
#' @param x An object returned by \code{summary}
#' @param ... Other parameters to be passed to \code{print.sumlucid_parallel}
#' @return \code{x}, invisibly. Called for its side effect: printing a structured
#'   parallel-model summary -- per-layer missing-data profile, overall and
#'   per-layer feature-selection overview, model fit statistics, regularization
#'   settings, and parameter estimates for sections (1) Y, (2) Z and (3) E. If
#'   \code{boot.se} was provided to \code{summary()}, bootstrap confidence limits
#'   are shown for each of those three sections.
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
      print(.append_significance(se_mu[[i]]))
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
  .e_src_par <- if (in_serial && length(transition_labels) > 0) {
    "each cluster from the previous serial stage (for each layer)"
  } else {
    "each exposure for each layer"
  }
  if (is.null(x$boot.se) || is.null(x$boot.se$beta)) {
    cat(sprintf("(%d) E: intercept and odds ratio of being assigned to each latent cluster for %s\n",
                e_idx, .e_src_par))
  } else {
    cat(sprintf("(%d) E: cluster-assignment coefficients (log-odds scale) for %s, with bootstrap interval and the odds ratio (exp of each column). `sig` marks intervals excluding the null.\n",
                e_idx, .e_src_par))
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

      # Any beta_use columns past the exposure block are CoG covariates -- always
      # shown (they are not penalty-selected), so the (3) E table is the same
      # with and without a bootstrap CI, and across model types.
      cov_names <- character(0)
      if (length(transition_labels) == 0 && (ncol(beta_use) - 1L) > n_exposure) {
        cov_names <- colnames(beta_use)[seq.int(n_exposure + 2L, ncol(beta_use))]
        g_names <- c(g_names, cov_names)
        selected_g <- c(selected_g, rep(TRUE, length(cov_names)))
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
            # keep the intercept, the selected exposures, and every covariate
            # row (a covariate name is not among the exposure names `g_names`)
            keep <- feature %in% selected_features | !(feature %in% g_names)
            beta_ci <- beta_ci[keep, , drop = FALSE]
          }
        }
        print(.with_or_columns(.append_significance(beta_ci)))
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
#' @return \code{x}, invisibly. Called for its side effect: printing the serial
#'   model's summary stage by stage -- per-stage missing-data profile,
#'   feature-selection overview, model fit statistics and parameter estimates,
#'   with the outcome section attached to the final stage. If \code{boot.se} was
#'   supplied to \code{summary()}, bootstrap confidence limits are shown
#'   alongside the estimates.
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



#' Print the outcome section of \code{print.sumlucid_early} for a normal Y
#'
#' @param x The \code{parameters$gamma} component of a
#'   \code{sumlucid_early} object.
#' @param K Number of clusters.
#' @param se Optional bootstrap CI data frame for gamma.
#' @return \code{x} invisibly printed as a side effect; called for its
#'   \code{print()} side effect.
#' @noRd
f.normal.early <- function(x, K, se){

  cat("(1) Y (continuous outcome): cluster 1 mean (intercept), mean differences for other clusters, and covariate effects \n")

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
  print(.append_significance(y))
}


#' Print the outcome section of \code{print.sumlucid_early} for a binary Y
#'
#' @param x The \code{parameters$gamma} component of a
#'   \code{sumlucid_early} object.
#' @param K Number of clusters.
#' @param se Optional bootstrap CI data frame for gamma.
#' @return Called for its \code{print()} side effect.
#' @noRd
f.binary.early <- function(x, K, se){
  if(is.null(se)){
    cat("(1) Y (binary outcome): log odds of Y for cluster 1 (reference) and log OR for rest cluster (and log OR of covariate if included)\n")
    gamma <- as.data.frame(x$beta)
    colnames(gamma) <- "gamma"
    gamma$`exp(gamma)` <- exp(gamma$gamma)
    print(.append_significance(gamma))
  } else {
    cat("(1) Y (binary outcome): log odds / log OR (coefficient scale) with bootstrap interval, and the odds ratio (exp of each column). `sig` marks intervals excluding the null.\n")
    # Consume the bootstrap CI table directly -- its `estimate` column already
    # is the point estimate -- instead of joining it to `x$beta` by row name.
    # A serial model's final stage prefixes the outcome coefficient names with
    # "Y." (extract_early_stage_vector), which no `x$beta` name would match, so
    # the join produced an all-NA table for serial + binary.
    gamma <- as.data.frame(se)
    rn <- rownames(gamma)
    rn <- sub("^Y\\.", "", rn)
    rn[tolower(rn) == "intercept"] <- "(Intercept)"
    rownames(gamma) <- rn
    if (ncol(gamma) >= 1 && identical(colnames(gamma)[1], "estimate")) {
      colnames(gamma)[1] <- "gamma"
    }
    gamma <- .append_significance(gamma)
    print(.with_or_columns(gamma))
  }
}

#' Append odds-ratio columns (exp of the coefficient-scale estimate and bounds)
#'
#' The bootstrap CI tables print on the coefficient (log-odds / linear-predictor)
#' scale so the null is a single 0 and \code{sig} is unambiguous. For a binary
#' outcome or a multinomial-logit assignment model those coefficients are log
#' odds ratios, so \code{exp()} of the estimate and of each interval bound is a
#' valid odds-ratio view; this adds it alongside without disturbing the
#' coefficient-scale columns or \code{sig}.
#'
#' @param df A data frame carrying a coefficient-scale \code{estimate}/first
#'   column and \code{norm_lower}/\code{norm_upper} bounds.
#' @return \code{df} with \code{OR}, \code{OR_lower}, \code{OR_upper} appended
#'   (before \code{sig} if present).
#' @noRd
.with_or_columns <- function(df) {
  df <- as.data.frame(df, check.names = FALSE)
  est_col <- if ("estimate" %in% names(df)) "estimate" else
    if ("gamma" %in% names(df)) "gamma" else if ("Gamma" %in% names(df)) "Gamma" else names(df)[1]
  or <- data.frame(OR = exp(df[[est_col]]), check.names = FALSE)
  if ("norm_lower" %in% names(df)) or$OR_lower <- exp(df[["norm_lower"]])
  if ("norm_upper" %in% names(df)) or$OR_upper <- exp(df[["norm_upper"]])
  if ("sig" %in% names(df)) {
    cbind(df[, setdiff(names(df), "sig"), drop = FALSE], or, sig = df[["sig"]])
  } else {
    cbind(df, or)
  }
}

#' Translate internal parallel-model outcome coefficient names for display
#'
#' Converts names like \code{"LC3"}/\code{"gamma2"} into
#' \code{"cluster3Layer1"}-style labels naming the layer and non-reference
#' cluster they refer to.
#'
#' @param nm Internal coefficient names, or \code{NULL} for generic
#'   fallback names.
#' @param K Per-layer cluster counts.
#' @return A character vector of display names, same length as \code{nm}.
#' @noRd
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

#' Print the outcome section of \code{print.sumlucid_parallel} for a normal Y
#'
#' @param x The \code{parameters$gamma} component of a
#'   \code{sumlucid_parallel} object.
#' @param K Per-layer cluster counts.
#' @param se Optional bootstrap CI data frame for gamma.
#' @return Called for its \code{print()} side effect.
#' @noRd
f.normal.parallel <- function(x, K, se){

  cat("(1) Y (continuous outcome): intercept, effects of each non-reference latent cluster for each layer of Y (and effect of covariates if included) \n")

  if(!is.null(se)){
    gamma <- as.data.frame(se)
    if(ncol(gamma) >= 1 && identical(colnames(gamma)[1], "estimate")) {
      colnames(gamma)[1] <- "Gamma"
    }
    rownames(gamma) <- normalize_parallel_y_names(rownames(gamma), K)
  } else {
    coef_vec <- parallel_delta_coef(x$Gamma)
    gamma <- data.frame(Gamma = as.numeric(coef_vec), check.names = FALSE)
    rownames(gamma) <- normalize_parallel_y_names(names(coef_vec), K)
  }
  print(.append_significance(gamma))
}


#' Print the outcome section of \code{print.sumlucid_parallel} for a binary Y
#'
#' @param x The \code{parameters$gamma} component of a
#'   \code{sumlucid_parallel} object.
#' @param K Per-layer cluster counts.
#' @param se Optional bootstrap CI data frame for gamma.
#' @return Called for its \code{print()} side effect.
#' @noRd
f.binary.parallel <- function(x, K, se){
  coef_vec <- parallel_delta_coef(x$Gamma)
  gamma <- data.frame(gamma = as.numeric(coef_vec), check.names = FALSE)
  rownames(gamma) <- normalize_parallel_y_names(names(coef_vec), K)
  if(is.null(se)){
    cat("(1) Y (binary outcome): intercept and log OR for non-reference clusters for each layer (and log OR of covariate if included)\n")
    gamma$`exp(gamma)` <- exp(gamma$gamma)
    print(.append_significance(gamma))
  } else {
    cat("(1) Y (binary outcome): log odds / log OR (coefficient scale) with bootstrap interval, and the odds ratio (exp of each column). `sig` marks intervals excluding the null.\n")
    se_df <- as.data.frame(se)
    rownames(se_df) <- normalize_parallel_y_names(rownames(se_df), K)
    idx <- match(rownames(gamma), rownames(se_df))
    se_df <- se_df[idx, , drop = FALSE]
    gamma <- cbind(gamma, se_df[, setdiff(colnames(se_df), "estimate"), drop = FALSE])
    print(.with_or_columns(.append_significance(gamma)))
  }
}

##########functions for LUCID in parallel##########

#' BIC for a parallel-model fit (Eqs 13/18, via \code{\link{parallel_npars}})
#'
#' @param object A fitted \code{lucid_parallel} object.
#' @return The scalar BIC.
#' @noRd
cal_bic_parallel <- function(object) {
  ##not having regularity yet, to be added
  s1 <- object$select$selectG
  s2 <- object$select$selectZ
  nG <- if (is.list(s1)) sum(sapply(s1, sum)) else sum(s1)
  nZ <- sapply(s2, function(x) sum(.selectZ_feature_mask(x)))
  npars <- parallel_npars(object)

  BIC <- -2 * object$likelihood + npars * log(nrow(object$inclusion.p[[1]]))

  return(BIC)
}

#' BIC for a serial-model fit: the sum of each submodel's own BIC
#'
#' A serial chain is a sequence of conditionally-fitted stages, not one
#' joint model, so there is no single joint BIC formula -- this sums each
#' stage's own BIC (via \code{\link{summary_lucid_simple}}), matching how
#' \code{\link{cal_loglik_serial}} sums log-likelihoods.
#'
#' @param object A fitted \code{lucid_serial} object.
#' @return The scalar sum of per-stage BICs.
#' @noRd
cal_bic_serial <- function(object) {

  sum_all_sub = summary_lucid_simple(object)
  BIC = 0
  for (i in 1:length(sum_all_sub)){
    BIC_temp = sum_all_sub[[i]]$BIC
    BIC = BIC + BIC_temp
    }
  return(BIC)
}

#' Log-likelihood for a serial-model fit: the sum of each submodel's own
#'
#' See \code{\link{cal_bic_serial}} for why this is a sum rather than a
#' single joint value. This is what \code{estimate_lucid()}'s serial branch
#' surfaces as the top-level \code{fit$likelihood}.
#'
#' @param object A fitted \code{lucid_serial} object.
#' @return The scalar sum of per-stage log-likelihoods.
#' @noRd
cal_loglik_serial <- function(object) {

  sum_all_sub = summary_lucid_simple(object)
  loglik = 0
  for (i in 1:length(sum_all_sub)){
    loglik_temp = sum_all_sub[[i]]$loglik
    loglik = loglik + loglik_temp
  }
  return(loglik)
}

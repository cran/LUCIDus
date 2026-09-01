# User-facing extractor functions that take a fitted LUCID model (any of the
# three model types, auto-detected from its class) and pull out one specific
# thing: selected exposures, selected omics features, hard cluster
# assignment, or the top-N most important omics features. None of these
# compute anything new -- they normalize fields/logic that already exist on
# the fitted object (or, for feature importance, in
# `plot_cluster_omic_profile.R`) into one consistent, class-generic entry
# point, so a user never has to branch on model type themselves.

#' Detect which LUCID model type a fitted object is
#'
#' @param model A fitted LUCID model.
#' @return "early", "parallel", or "serial".
#' @noRd
.detect_lucid_model <- function(model) {
  if (inherits(model, "early_lucid")) return("early")
  if (inherits(model, "lucid_parallel")) return("parallel")
  if (inherits(model, "lucid_serial")) return("serial")
  stop("`model` must be a fitted LUCID model (early_lucid, lucid_parallel, ",
      "or lucid_serial), as returned by estimate_lucid() or lucid(); got ",
      "class '", paste(class(model), collapse = "/"), "'.", call. = FALSE)
}

#' Collapse a selectZ entry to one logical value per feature
#'
#' A layer's \code{selectZ} is either a plain logical vector (one entry per
#' feature) or a \code{K x M} logical matrix (one entry per cluster per
#' feature, for a per-cluster selection); this collapses either shape to one
#' logical value per feature, \code{TRUE} if the feature is selected in at
#' least one cluster.
#'
#' @param x A logical vector, a logical matrix, or \code{NULL}.
#' @return A logical vector (one entry per feature), or \code{NULL} if
#'   \code{x} is \code{NULL}.
#' @noRd
.selectZ_feature_mask <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.null(dim(x))) as.logical(x) else colSums(x) > 0
}

#' Recursively compute the hard (argmax) cluster assignment
#'
#' Walks the same matrix-or-list structure as \code{inclusion.p}: a matrix
#' leaf becomes one hard assignment per row (via \code{nnet::which.is.max()},
#' matching \code{predict_lucid()}'s own tie-breaking convention exactly); a
#' list recurses.
#'
#' @param p An \code{inclusion.p}-shaped matrix or (possibly nested) list of
#'   matrices.
#' @return The same shape as \code{p}, with each matrix leaf replaced by a
#'   numeric vector of hard cluster labels (\code{1, ..., K}).
#' @noRd
.hard_assign_recursive <- function(p) {
  if (is.list(p)) {
    return(lapply(p, .hard_assign_recursive))
  }
  p <- as.matrix(p)
  as.numeric(vapply(seq_len(nrow(p)), function(i) nnet::which.is.max(p[i, ]),
                    numeric(1)))
}

#' Extract selected (retained) exposures from a fitted LUCID model
#'
#' Reads \code{model$select$selectG} (or, for a serial model, stage 1's own
#' selection), auto-detecting the model type from \code{class(model)}. Stage
#' 1 is the only stage in a serial chain whose "G" is the cohort's actual
#' exposures -- from stage 2 on, "G" is the previous stage's posterior
#' cluster probabilities, so there is nothing there for an exposure-selection
#' result to describe (see \code{\link{estimate_lucid}}'s \code{@return} for
#' the full rationale). This is exactly the same field \code{estimate_lucid()}
#' already returns; this function only adds the class-based dispatch and the
#' union/per-layer/per-stage bookkeeping so a caller doesn't have to.
#'
#' @param model A fitted \code{early_lucid}, \code{lucid_parallel}, or
#'   \code{lucid_serial} object.
#' @param layer For a parallel model only: which layer's own exposure
#'   selection to return (an integer index or a layer name). If \code{NULL}
#'   (the default), returns the union across layers (selected in at least one
#'   layer) -- ignored for early and serial.
#' @return A named logical vector, one entry per exposure, \code{TRUE} where
#'   retained.
#' @examples
#' idx <- 1:200
#' G <- sim_data$G[idx, ]
#' Z <- sim_data$Z[idx, ]
#' Y_normal <- sim_data$Y_normal[idx, ]
#' fit <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early",
#'                       family = "normal", K = 2, Rho_G = 0.1,
#'                       max_itr = 10, max_tot.itr = 30)
#' get_selected_G(fit)
#' @export
get_selected_G <- function(model, layer = NULL) {
  lucid_model <- .detect_lucid_model(model)
  if (lucid_model == "serial") {
    # stage 1's own selection only -- see the roxygen above for why
    model <- model$submodel[[1]]
    lucid_model <- .detect_lucid_model(model)
  }
  mask <- if (lucid_model == "parallel" && !is.null(layer)) {
    model$select$selectG_layer[[layer]]
  } else {
    model$select$selectG
  }
  # var.names$Gnames is the exposures followed by any CoG covariate names
  # (estimate_lucid() builds it as `c(Gnames, CoGnames)`); selectG only ever
  # covers the exposure portion, so take just the leading names that match.
  stats::setNames(mask, model$var.names$Gnames[seq_along(mask)])
}

#' Extract selected (retained) omics features from a fitted LUCID model
#'
#' Reads \code{model$select$selectZ}, auto-detecting the model type from
#' \code{class(model)} and collapsing any per-cluster selection matrix (a
#' parallel-model layer's \code{selectZ} can be a \code{K x M} matrix rather
#' than a plain vector) to one logical value per feature via an internal
#' helper. Unlike exposure selection, every
#' stage of a serial model has a meaningful omics selection, so this returns
#' a per-stage breakdown rather than one stage's alone.
#'
#' @param model A fitted \code{early_lucid}, \code{lucid_parallel}, or
#'   \code{lucid_serial} object.
#' @param layer For a parallel model (or a serial stage that is itself a
#'   parallel sub-model): which layer's own omics selection to return. If
#'   \code{NULL} (the default), returns a named list, one entry per layer.
#' @param stage For a serial model only: which stage's own omics selection to
#'   return. If \code{NULL} (the default), returns a named list, one entry
#'   per stage (each shaped like this function's early/parallel return,
#'   depending on that stage's own type).
#' @return A named logical vector (early; parallel with \code{layer} given),
#'   a named list of logical vectors (parallel with \code{layer = NULL}), or
#'   a named list of per-stage results (serial).
#' @examples
#' idx <- 1:200
#' G <- sim_data$G[idx, ]
#' Z <- sim_data$Z[idx, ]
#' Y_normal <- sim_data$Y_normal[idx, ]
#' fit <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early",
#'                       family = "normal", K = 2, Rho_Z_Mu = 5,
#'                       max_itr = 10, max_tot.itr = 30)
#' get_selected_Z(fit)
#' @export
get_selected_Z <- function(model, layer = NULL, stage = NULL) {
  lucid_model <- .detect_lucid_model(model)

  extract_early_or_parallel <- function(m, layer) {
    m_type <- .detect_lucid_model(m)
    if (m_type == "early") {
      return(stats::setNames(.selectZ_feature_mask(m$select$selectZ),
                             m$var.names$Znames))
    }
    # parallel
    if (!is.null(layer)) {
      return(stats::setNames(.selectZ_feature_mask(m$select$selectZ[[layer]]),
                             m$var.names$Znames[[layer]]))
    }
    masks <- lapply(seq_along(m$select$selectZ), function(i) {
      stats::setNames(.selectZ_feature_mask(m$select$selectZ[[i]]),
                      m$var.names$Znames[[i]])
    })
    names(masks) <- .layer_labels(m$Z, length(masks))
    masks
  }

  if (lucid_model != "serial") {
    return(extract_early_or_parallel(model, layer))
  }

  # serial: per-stage breakdown, each shaped like the early/parallel result
  submodels <- model$submodel
  out <- lapply(submodels, extract_early_or_parallel, layer = layer)
  names(out) <- .layer_labels(model$Z, length(out), "Stage")
  if (!is.null(stage)) return(out[[stage]])
  out
}

#' Extract the hard cluster assignment from a fitted LUCID model
#'
#' Computes the maximum-a-posteriori cluster label for every observation
#' directly from \code{model$inclusion.p}, with no need to re-run prediction
#' or supply \code{G}/\code{Z}/\code{Y} again. Shaped exactly like
#' \code{predict_lucid()}'s own \code{pred.x}: a numeric vector for early, a
#' list by layer for parallel, and a list by stage (recursively shaped) for
#' serial. Labels run \code{1, ..., K}, matching Eq 21 and the row names used
#' by \code{summary()}.
#'
#' @param model A fitted \code{early_lucid}, \code{lucid_parallel}, or
#'   \code{lucid_serial} object.
#' @return A numeric vector (early), a named list by layer (parallel), or a
#'   named list by stage (serial, each element shaped like the above
#'   depending on that stage's own type).
#' @examples
#' idx <- 1:200
#' G <- sim_data$G[idx, ]
#' Z <- sim_data$Z[idx, ]
#' Y_normal <- sim_data$Y_normal[idx, ]
#' fit <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early",
#'                       family = "normal", K = 2,
#'                       max_itr = 10, max_tot.itr = 30)
#' table(get_cluster_assignment(fit))
#' @export
get_cluster_assignment <- function(model) {
  .detect_lucid_model(model)  # validates model's class; result unused here
  .hard_assign_recursive(model$inclusion.p)
}

#' Extract the top-N most important omics features from a fitted LUCID model
#'
#' Reuses \code{plot_cluster_omic_profile()}'s own feature-ranking criterion
#' (see its documentation for what "separation" means) rather than
#' introducing a second ranking rule: this is the same score
#' \code{plot_cluster_omic_profile()} sorts features by, just returned as
#' data instead of a plot. One panel is produced per relevant unit -- the
#' whole omics matrix for early, one per layer for parallel, and for serial
#' one per stage (or one per layer within a stage that is itself parallel).
#'
#' @param model A fitted \code{early_lucid}, \code{lucid_parallel}, or
#'   \code{lucid_serial} object.
#' @param top_n Number of top features to return per panel (default 10). If
#'   a panel has fewer features than \code{top_n}, all of them are returned.
#' @param importance Ranking criterion: "separation" (between-cluster spread
#'   over within-cluster SD, the default), "range", or "sd" of the cluster
#'   means -- identical meaning to \code{plot_cluster_omic_profile()}'s own
#'   \code{importance} argument.
#' @return A named list, one entry per panel (layer/stage), each a named
#'   numeric vector of the top \code{top_n} features by \code{importance},
#'   sorted descending.
#' @examples
#' idx <- 1:200
#' G <- sim_data$G[idx, ]
#' Z <- sim_data$Z[idx, ]
#' Y_normal <- sim_data$Y_normal[idx, ]
#' fit <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early",
#'                       family = "normal", K = 2,
#'                       max_itr = 10, max_tot.itr = 30)
#' get_top_omics_features(fit, top_n = 3)
#' @export
get_top_omics_features <- function(model, top_n = 10,
                                   importance = c("separation", "range", "sd")) {
  .detect_lucid_model(model)  # validates model's class
  importance <- match.arg(importance)
  if (!is.numeric(top_n) || length(top_n) != 1 || is.na(top_n) || top_n < 1) {
    stop("`top_n` must be a single positive integer.", call. = FALSE)
  }
  panels <- .lucid_omic_panels(model)
  lapply(panels, function(p) {
    score <- .feature_importance(p, importance)
    n_keep <- min(as.integer(top_n), length(score))
    sort(score, decreasing = TRUE)[seq_len(n_keep)]
  })
}

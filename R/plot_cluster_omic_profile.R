# Per-cluster omics profiles for early, parallel and serial fits.
#
# LUCID's clusters are defined by their omics profiles, but until now nothing in
# the package showed those profiles: plot() draws the path structure, and
# summary() prints res_Mu as a numeric table. This file answers the question a
# reader actually has -- which features distinguish the clusters, and how.
#
# The three model types store the same information in three different shapes, so
# everything here works on "panels" produced by .lucid_omic_panels(): once a fit
# has been reduced to a list of (mu, var, label), the ranking and both renderers
# are model-agnostic.

# ggplot2 refers to data columns by bare name; declare them so R CMD check
# does not report them as undefined globals.
utils::globalVariables(c("cluster", "feature", "value"))

# ---- panel extraction -------------------------------------------------------

#' Within-cluster variance per (cluster, feature), from either model's shape
#'
#' Early keeps a list of \code{K} covariance matrices; parallel keeps an
#' \code{M x M x K} array. \code{simplify2array()} collapses the first onto
#' the second. Falls back to unit variance for an unrecognized shape, so a
#' profile can still be drawn with the unstandardized measures.
#'
#' @param sigma Either shape of covariance storage.
#' @param K Number of clusters.
#' @param M Number of features.
#' @return A \code{K x M} matrix of within-cluster variances.
#' @noRd
.panel_var <- function(sigma, K, M) {
  arr <- if (is.list(sigma)) simplify2array(sigma) else sigma
  if (!is.array(arr) || length(dim(arr)) != 3L) {
    # A shape we do not recognise: fall back to unit variance rather than fail,
    # so a profile can still be drawn with the unstandardized measures.
    return(matrix(1, nrow = K, ncol = M))
  }
  out <- matrix(NA_real_, nrow = K, ncol = M)
  for (k in seq_len(K)) out[k, ] <- diag(arr[, , k])
  out
}

#' Build one plotting panel from a mu matrix plus its covariance object
#'
#' @param mu Cluster means (K x M matrix).
#' @param sigma Cluster covariances, either shape \code{\link{.panel_var}}
#'   accepts.
#' @param label This panel's display label.
#' @return A list: \code{mu}, \code{var} (from \code{.panel_var()}, with
#'   \code{mu}'s feature names carried over), \code{label}.
#' @noRd
.make_panel <- function(mu, sigma, label) {
  mu <- as.matrix(mu)
  K <- nrow(mu); M <- ncol(mu)
  if (is.null(colnames(mu))) colnames(mu) <- paste0("feature", seq_len(M))
  v <- .panel_var(sigma, K, M)
  # carry mu's feature names onto the variance matrix so both can be indexed by
  # feature name once the top-n subset is chosen
  dimnames(v) <- list(NULL, colnames(mu))
  list(mu = mu, var = v, label = label)
}

#' Recover layer/stage names, falling back to numbered labels
#'
#' Layer names survive only on \code{fit$Z}, never on \code{res_Mu} or
#' \code{var.names$Znames}.
#'
#' @param z The fit's \code{Z} (or a stage's), whose names (if any) are used.
#' @param n Number of layers/stages, for the fallback numbering.
#' @param prefix Fallback label prefix (default \code{"Layer"}).
#' @return A character vector of length \code{n}.
#' @noRd
.layer_labels <- function(z, n, prefix = "Layer") {
  nm <- names(z)
  if (is.null(nm) || any(!nzchar(nm))) paste(prefix, seq_len(n)) else nm
}

#' Reduce any fitted LUCID model to a flat list of omics panels
#'
#' Each panel is one cluster-by-feature block that can be plotted on its own:
#' the whole omics matrix for an early fit, one per layer for a parallel fit,
#' and for a serial fit one per stage -- or one per layer within a stage that is
#' itself a parallel sub-model.
#'
#' @param x a fitted LUCID model
#' @return a named list of panels, each with `mu`, `var` and `label`
#' @noRd
.lucid_omic_panels <- function(x) {
  if (inherits(x, "early_lucid")) {
    lab <- .layer_labels(x$Z, 1L, "Omics")
    if (identical(lab[1], "Omics 1")) lab[1] <- "Omics"
    p <- list(.make_panel(x$res_Mu, x$res_Sigma, lab[1]))
    names(p) <- lab[1]
    return(p)
  }

  if (inherits(x, "lucid_parallel")) {
    n <- length(x$res_Mu)
    labs <- .layer_labels(x$Z, n)
    p <- lapply(seq_len(n), function(i) {
      .make_panel(x$res_Mu[[i]], x$res_Sigma[[i]], labs[i])
    })
    names(p) <- labs
    return(p)
  }

  if (inherits(x, "lucid_serial")) {
    stage_labs <- .layer_labels(x$Z, length(x$submodel), "Stage")
    out <- list()
    for (i in seq_along(x$submodel)) {
      sm <- x$submodel[[i]]
      if (inherits(sm, "lucid_parallel")) {
        # a stage that is itself a parallel sub-model contributes one panel per
        # layer, labelled with both so the reader can place it in the chain
        inner <- .layer_labels(sm$Z, length(sm$res_Mu))
        for (j in seq_along(sm$res_Mu)) {
          lab <- paste0(stage_labs[i], " - ", inner[j])
          out[[lab]] <- .make_panel(sm$res_Mu[[j]], sm$res_Sigma[[j]], lab)
        }
      } else {
        lab <- stage_labs[i]
        out[[lab]] <- .make_panel(sm$res_Mu, sm$res_Sigma, lab)
      }
    }
    return(out)
  }

  stop("`x` must be a fitted LUCID model (early_lucid, lucid_parallel or ",
       "lucid_serial).", call. = FALSE)
}

# ---- feature importance -----------------------------------------------------

#' Rank features within one panel by importance
#'
#' See the roxygen on \code{plot_cluster_omic_profile()} for why
#' "separation" is the default.
#'
#' @param panel A panel from \code{\link{.make_panel}}.
#' @param importance "separation" (between-cluster spread over within-cluster
#'   SD), "range", or "sd" (of the cluster means).
#' @return A named numeric vector of importance scores, one per feature.
#' @noRd
.feature_importance <- function(panel, importance = c("separation", "range", "sd")) {
  importance <- match.arg(importance)
  mu <- panel$mu
  spread <- switch(importance,
    separation = apply(mu, 2, stats::sd),
    sd         = apply(mu, 2, stats::sd),
    range      = apply(mu, 2, function(v) diff(range(v)))
  )
  if (importance == "separation") {
    within <- sqrt(pmax(colMeans(panel$var), .Machine$double.eps))
    spread <- spread / within
  }
  stats::setNames(as.numeric(spread), colnames(mu))
}

# ---- shared data assembly ---------------------------------------------------

#' Assemble one panel's plotting data frame (top features, scaled, labeled)
#'
#' Selects the top \code{top_n} features by importance, optionally z-scores
#' (or, at exactly \code{K = 2} with the default \code{scale}, centres
#' without rescaling -- a z-score over two values is always exactly
#' \eqn{\pm 1/\sqrt{2}}, which would saturate every tile and show sign but
#' no magnitude), and reshapes to one row per (feature, cluster).
#'
#' @param panel A panel from \code{\link{.make_panel}}.
#' @param top_n Number of features to keep (fewer available is not an error).
#' @param importance Passed to \code{\link{.feature_importance}}.
#' @param scale Whether to z-score each feature across clusters.
#' @param cluster_labels Cluster display labels, one per cluster.
#' @param scale_is_default Whether \code{scale} is still at its default
#'   (controls the K = 2 centre-only special case).
#' @return A data frame: \code{feature}, \code{cluster}, \code{value}
#'   (plotted fill/bar value), \code{mean}, \code{sd}, \code{score}, with
#'   a \code{"centred"} attribute recording whether the K = 2 case applied.
#' @noRd
.panel_frame <- function(panel, top_n, importance, scale, cluster_labels,
                         scale_is_default = FALSE) {
  score <- .feature_importance(panel, importance)
  # fewer features than requested is not an error, it just means show them all
  n_keep <- min(top_n, length(score))
  keep <- names(sort(score, decreasing = TRUE))[seq_len(n_keep)]

  mu <- panel$mu[, keep, drop = FALSE]
  K <- nrow(mu)
  clusters <- if (is.null(cluster_labels)) paste("Cluster", seq_len(K)) else cluster_labels
  if (length(clusters) != K) {
    stop("`cluster_labels` must have one entry per cluster (", K, ").", call. = FALSE)
  }

  # z-score each feature across clusters: the single-cell convention, and what
  # makes features on different scales comparable inside one panel.
  #
  # With only two clusters that standardization is degenerate -- a z-score over
  # two values is always exactly +/-1/sqrt(2) -- so every tile would saturate and
  # the fill would carry sign but no magnitude. At K = 2 the default therefore
  # centres without rescaling, which keeps the palette diverging while letting
  # the size of each difference show. An explicit scale = TRUE is still honoured.
  centre_only <- isTRUE(scale) && nrow(mu) == 2L && isTRUE(scale_is_default)
  fill <- mu
  if (centre_only) {
    fill <- sweep(mu, 2, colMeans(mu), "-")
  } else if (isTRUE(scale)) {
    fill <- apply(mu, 2, function(v) {
      s <- stats::sd(v)
      if (!is.finite(s) || s < .Machine$double.eps) rep(0, length(v)) else (v - mean(v)) / s
    })
    fill <- matrix(fill, nrow = K, dimnames = dimnames(mu))
  }

  out <- data.frame(
    feature = factor(rep(keep, each = K), levels = rev(keep)),
    cluster = factor(rep(clusters, times = length(keep)), levels = clusters),
    value   = as.numeric(fill),
    mean    = as.numeric(mu),
    sd      = as.numeric(sqrt(panel$var[, keep, drop = FALSE])),
    score   = rep(score[keep], each = K),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(out, "centred") <- centre_only
  out
}

#' Shared ggplot2 theme for the cluster-omics-profile plots
#'
#' @return A \code{ggplot2} theme object.
#' @noRd
.profile_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(size = 9),
      axis.text.x      = ggplot2::element_text(size = 9),
      plot.title       = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey30"),
      plot.caption     = ggplot2::element_text(size = 8, colour = "grey45"),
      legend.position  = "right"
    )
}

#' Default per-layer hues, cycled if there are more layers than colours
#'
#' The first two match the package's own omics and cluster colours from
#' \code{plot_lucid.R}.
#'
#' @param n Number of layers.
#' @return A character vector of \code{n} hex colours.
#' @noRd
.default_layer_colors <- function(n) {
  base <- c("#2fa4da", "#eb8c30", "#67928b", "#9d7bb0", "#c4622d", "#4f7871")
  rep(base, length.out = n)
}

# The heatmap fill is signed, so it uses a fixed diverging scale rather than the
# layer hue: cool for low, warm for high, which is the convention a reader of
# single-cell marker heatmaps already knows. Deriving the scale from the layer
# hue instead would invert that convention for any layer whose hue is warm, so
# the layer is identified by its subtitle colour and by the bar plot instead.
.HEAT_LOW  <- "#4575b4"
.HEAT_HIGH <- "#d73027"

#' Legend title for the plotted value, matching what \code{.panel_frame} computed
#'
#' @param scale Whether the value is a z-score (vs. a raw cluster mean).
#' @param wrap Whether to insert a line break for a narrower legend.
#' @param centred Whether the K = 2 centre-only case applied (see
#'   \code{\link{.panel_frame}}).
#' @return A character string.
#' @noRd
.value_label <- function(scale, wrap = TRUE, centred = FALSE) {
  if (!isTRUE(scale)) return("cluster mean")
  if (centred) return(if (wrap) "difference from\nfeature mean" else "difference from feature mean")
  if (wrap) "z-score\nacross clusters" else "z-score across clusters"
}

#' Plot caption describing the feature-ranking criterion
#'
#' @param importance As passed to \code{\link{.feature_importance}}.
#' @return A character string.
#' @noRd
.importance_caption <- function(importance) {
  switch(importance,
    separation = "Features ranked by between-cluster spread / within-cluster SD",
    range      = "Features ranked by the range of cluster means",
    sd         = "Features ranked by the SD of cluster means"
  )
}

#' Mix a colour with white
#'
#' \code{grey85} as a ramp start leaves the palest cluster nearly invisible
#' against the panel; a tint of the layer hue stays legible and keeps the
#' whole ramp within that layer's identity.
#'
#' @param hue A colour, as accepted by \code{grDevices::col2rgb()}.
#' @param amount Fraction of the way to white (default \code{0.78}).
#' @return A hex colour string.
#' @noRd
.tint <- function(hue, amount = 0.78) {
  rgb <- grDevices::col2rgb(hue)[, 1]
  mixed <- rgb + (255 - rgb) * amount
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

# ---- renderers --------------------------------------------------------------

#' Render one panel as a heatmap
#'
#' @param df A panel data frame from \code{\link{.panel_frame}}.
#' @param label This panel's subtitle.
#' @param importance Passed to \code{\link{.importance_caption}}.
#' @param scale Passed to \code{\link{.value_label}}.
#' @param hue This panel's layer colour (used for the subtitle).
#' @return A \code{ggplot} object.
#' @noRd
.profile_heatmap <- function(df, label, importance, scale, hue) {
  n_clust <- nlevels(df$cluster)
  lim <- max(abs(df$value), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = feature, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(
      low = .HEAT_LOW, mid = "white", high = .HEAT_HIGH, midpoint = 0,
      limits = c(-lim, lim), name = .value_label(scale, centred = isTRUE(attr(df, "centred")))
    ) +
    ggplot2::labs(
      title    = "Cluster omics profile",
      subtitle = label,
      caption  = .importance_caption(importance),
      x = NULL, y = NULL
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    .profile_theme() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 10, colour = hue, face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = if (n_clust > 6) 45 else 0,
        hjust = if (n_clust > 6) 0 else 0.5
      )
    )
}

#' Render one panel as a bar plot
#'
#' @param df A panel data frame from \code{\link{.panel_frame}}.
#' @param label This panel's subtitle.
#' @param importance Passed to \code{\link{.importance_caption}}.
#' @param scale Passed to \code{\link{.value_label}}.
#' @param hue This panel's layer colour, ramped by cluster (via
#'   \code{\link{.tint}}).
#' @return A \code{ggplot} object.
#' @noRd
.profile_bar <- function(df, label, importance, scale, hue) {
  n_clust <- nlevels(df$cluster)
  # a sequential ramp within the layer's hue, so any number of clusters works --
  # a categorical palette would run out and start recycling around 8-12
  shades <- grDevices::colorRampPalette(c(.tint(hue), hue))(n_clust)

  ggplot2::ggplot(df, ggplot2::aes(x = value, y = feature, fill = cluster)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = shades, name = NULL) +
    ggplot2::labs(
      title    = "Cluster omics profile",
      subtitle = label,
      caption  = .importance_caption(importance),
      x = .value_label(scale, wrap = FALSE, centred = isTRUE(attr(df, "centred"))), y = NULL
    ) +
    .profile_theme() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 10, colour = hue, face = "bold")
    )
}

# ---- exported ---------------------------------------------------------------

#' @title Plot per-cluster omics profiles
#' @description
#' Shows which omics features distinguish the latent clusters, and in which
#' direction, for an early, parallel or serial fit. A parallel or serial model
#' produces one plot per omics layer, returned as a named list, so no figure has
#' to accommodate every layer at once.
#'
#' The default rendering is a cluster-by-feature heatmap in the style used for
#' single-cell cluster markers: features on the vertical axis ordered by how
#' strongly they separate the clusters, clusters across the top, and fill
#' showing how high or low each cluster sits for that feature.
#'
#' @section Which features are shown:
#' Only the `top_n` most discriminating features per panel are drawn, ranked by
#' `importance`:
#'
#' \describe{
#'   \item{`"separation"` (default)}{The spread of the cluster means divided by
#'     the typical within-cluster spread,
#'     \eqn{\mathrm{sd}_k(\mu_{kj}) / \sqrt{\overline{\sigma^2_{kj}}}}. This is
#'     the only option that accounts for noise: a feature whose cluster means
#'     differ by two units is uninformative if its within-cluster standard
#'     deviation is also two. It is scale-free, so features on different scales
#'     are comparable, and it is defined the same way for any number of
#'     clusters.}
#'   \item{`"range"`}{\eqn{\max_k \mu_{kj} - \min_k \mu_{kj}}, in the data's own
#'     units. The most directly interpretable, and the quantity the package's
#'     own feature selection thresholds.}
#'   \item{`"sd"`}{The standard deviation of the cluster means, without
#'     standardizing by within-cluster spread.}
#' }
#'
#' A feature the model deselected has identical means across clusters and so
#' scores zero under all three, sorting last.
#'
#' @section What the colour means:
#' With `scale = TRUE` (the default) the fill is a z-score computed across
#' clusters within each feature, so the palette is centred at zero and a
#' feature's own baseline does not dominate. With `scale = FALSE` the fill is
#' the fitted cluster mean on the data's original scale.
#'
#' **Two clusters are a special case.** A z-score over two values is always
#' exactly \eqn{\pm 1/\sqrt{2}}, so at `K = 2` every tile would saturate and
#' the fill would carry the sign of the difference but nothing about its size.
#' When `K = 2` and `scale` is left at its default, the fill is therefore
#' centred without rescaling -- each cluster mean minus that feature's average
#' across clusters -- which keeps the palette diverging while restoring
#' magnitude. Passing `scale = TRUE` explicitly overrides this and gives the
#' degenerate z-score. The legend always names whichever quantity was used, and
#' the returned data records it in the `"centred"` attribute.
#'
#' The underlying means and within-cluster standard deviations are always
#' available on the returned object, see *Value*.
#'
#' @param x A fitted LUCID model: `early_lucid`, `lucid_parallel` or
#'   `lucid_serial`.
#' @param type `"heatmap"` (default) or `"bar"`.
#' @param top_n Number of features to show per panel, default 10. A layer with
#'   fewer features than this shows all of them rather than erroring.
#' @param importance Ranking measure; see *Which features are shown*.
#' @param layer_names Character vector of layer or stage names, used as plot
#'   subtitles. Defaults to the names of the omics list the model was fitted to,
#'   falling back to `"Layer 1"`, `"Layer 2"` and so on.
#' @param layer_colors One colour per layer. The heatmap uses it as the high end
#'   of its diverging scale, and the bar plot as the darkest of a sequential
#'   ramp across clusters -- which is what lets the bar plot handle any number
#'   of clusters without recycling colours.
#' @param cluster_labels Optional labels for the clusters, defaulting to
#'   `"Cluster 1"`, `"Cluster 2"` and so on. Must have one entry per cluster.
#' @param scale If `TRUE` (default) the fill is a per-feature z-score across
#'   clusters; if `FALSE` it is the fitted cluster mean. At `K = 2` the default
#'   instead centres without rescaling, because a two-value z-score is
#'   degenerate -- see *What the colour means*.
#'
#' @return A named list of `ggplot` objects, one per omics layer -- length one
#'   for an early fit. Each carries the data it drew as the attribute
#'   `"profile_data"`: a data frame of `feature`, `cluster`, `value` (what is
#'   plotted), `mean` (the fitted cluster mean), `sd` (within-cluster standard
#'   deviation) and `score` (the importance value the ranking used), so the
#'   ranking can be extracted without recomputing it. That data frame also
#'   carries a `"centred"` attribute recording whether the fill was centred
#'   rather than z-scored, which is what happens by default at `K = 2`.
#'
#' @examples
#' \donttest{
#' # a small subset keeps the example quick
#' G <- sim_data$G[1:150, , drop = FALSE]
#' Z <- sim_data$Z[1:150, , drop = FALSE]
#' Y <- sim_data$Y_normal[1:150]
#'
#' fit <- estimate_lucid(G = G, Z = Z, Y = Y, lucid_model = "early",
#'                       family = "normal", K = 2, seed = 1008,
#'                       max_itr = 20, max_tot.itr = 50)
#'
#' p <- plot_cluster_omic_profile(fit)
#' p[[1]]
#'
#' # bar rendering, and more features
#' plot_cluster_omic_profile(fit, type = "bar", top_n = 15)[[1]]
#'
#' # the ranking behind the figure
#' head(unique(attr(p[[1]], "profile_data")[, c("feature", "score")]))
#' }
#'
#' @importFrom grDevices colorRampPalette
#' @export
plot_cluster_omic_profile <- function(x,
                                      type = c("heatmap", "bar"),
                                      top_n = 10,
                                      importance = c("separation", "range", "sd"),
                                      layer_names = NULL,
                                      layer_colors = NULL,
                                      cluster_labels = NULL,
                                      scale = TRUE) {
  type <- match.arg(type)
  importance <- match.arg(importance)
  if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1) {
    stop("`top_n` must be a single positive number.", call. = FALSE)
  }

  # whether the caller left `scale` alone decides the K = 2 fallback below
  scale_is_default <- missing(scale)

  panels <- .lucid_omic_panels(x)
  n <- length(panels)

  if (!is.null(layer_names)) {
    if (length(layer_names) != n) {
      stop("`layer_names` must have one entry per omics layer (", n, ").",
           call. = FALSE)
    }
    for (i in seq_len(n)) panels[[i]]$label <- layer_names[i]
    names(panels) <- layer_names
  }

  hues <- if (is.null(layer_colors)) .default_layer_colors(n) else
    rep(layer_colors, length.out = n)

  out <- lapply(seq_len(n), function(i) {
    df <- .panel_frame(panels[[i]], top_n, importance, scale, cluster_labels,
                       scale_is_default = scale_is_default)
    p <- if (type == "heatmap") {
      .profile_heatmap(df, panels[[i]]$label, importance, scale, hues[i])
    } else {
      .profile_bar(df, panels[[i]]$label, importance, scale, hues[i])
    }
    attr(p, "profile_data") <- df
    p
  })
  names(out) <- names(panels)
  out
}

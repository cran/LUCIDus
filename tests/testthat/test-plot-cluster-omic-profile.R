# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# plot_cluster_omic_profile(): the ranking, and the plumbing across model types.
#
# The ranking is tested on hand-built panels where mu and the within-cluster
# variances are known exactly, so the assertions are about the measure rather
# than about whatever a fit happened to produce. The plumbing is then tested on
# real fits, because the three model types store mu and Sigma in three different
# shapes and that is where this function is most likely to break.

# ---- the importance measure ------------------------------------------------

mk_panel <- function(mu, var, label = "test") {
  mu <- as.matrix(mu)
  if (is.null(colnames(mu))) colnames(mu) <- paste0("f", seq_len(ncol(mu)))
  v <- matrix(var, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  dimnames(v) <- list(NULL, colnames(mu))
  list(mu = mu, var = v, label = label)
}

test_that("separation ranks a well-separated feature above a noisy one", {
  # f1 and f2 have the SAME spread of cluster means, but f2 is four times as
  # noisy within cluster. Only a standardized measure can tell them apart -- and
  # telling them apart is the entire reason "separation" is the default.
  p <- mk_panel(mu = rbind(c(0, 0), c(2, 2)), var = c(0.25, 4))

  sep <- .feature_importance(p, "separation")
  rng <- .feature_importance(p, "range")

  expect_gt(sep[["f1"]], sep[["f2"]])
  expect_equal(unname(rng[["f1"]]), unname(rng[["f2"]]))   # range cannot
})

test_that("range and sd report the spread of the cluster means", {
  p <- mk_panel(mu = rbind(c(0, 0), c(3, 1)), var = c(1, 1))
  expect_equal(unname(.feature_importance(p, "range")), c(3, 1))
  expect_equal(unname(.feature_importance(p, "sd")),
               c(stats::sd(c(0, 3)), stats::sd(c(0, 1))))
})

test_that("a feature with identical cluster means scores zero", {
  # this is what a deselected feature looks like, so it must sort last
  p <- mk_panel(mu = rbind(c(1, 5), c(3, 5)), var = c(1, 1))
  for (m in c("separation", "range", "sd")) {
    s <- .feature_importance(p, m)
    expect_equal(unname(s[["f2"]]), 0, info = m)
    expect_gt(s[["f1"]], 0)
  }
})

test_that("zero within-cluster variance does not produce Inf or NaN", {
  p <- mk_panel(mu = rbind(c(0, 0), c(2, 2)), var = c(0, 1))
  s <- .feature_importance(p, "separation")
  expect_true(all(is.finite(s)))
})

# ---- panels across model types ---------------------------------------------

fit_of <- function(model, n_layer = 2L, K = 2, M = 8, seed = 21) {
  d <- sim_lucid(model, N = 250, K = K, M = M, n_layer = n_layer,
                 sep = 2.5, family = "normal", seed = seed)
  KK <- if (model == "early") K else as.list(rep(K, n_layer))
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                 family = "normal", K = KK, init_omic.data.model = "VVV")
  ))))
  fit
}

test_that("panels are extracted for every model type with the right count", {
  expect_length(.lucid_omic_panels(fit_of("early", 1L)), 1L)
  expect_length(.lucid_omic_panels(fit_of("parallel", 3L)), 3L)
  expect_length(.lucid_omic_panels(fit_of("serial", 2L)), 2L)
})

test_that("a serial stage that is itself parallel contributes one panel per layer", {
  topo <- list(2, list(2, 2))
  d <- sim_lucid_serial_topology(topo, N = 250, M = 6, sep = 2.5, seed = 22)
  names(d$Z) <- c("first", "second")
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                 family = "normal", K = topo, init_omic.data.model = "VVV")
  ))))
  p <- .lucid_omic_panels(fit)

  # stage 1 is early (one panel), stage 2 is parallel over two layers (two)
  expect_length(p, 3L)
  expect_match(names(p)[1], "first")
  expect_true(all(grepl("second", names(p)[2:3])))
})

test_that("layer names come from the omics list the model was fitted to", {
  d <- sim_lucid("parallel", N = 250, K = 2, M = 6, n_layer = 2, sep = 2.5, seed = 23)
  names(d$Z) <- c("methylome", "transcriptome")
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "parallel",
                 family = "normal", K = list(2, 2), init_omic.data.model = "VVV")
  ))))
  expect_equal(names(plot_cluster_omic_profile(fit)),
               c("methylome", "transcriptome"))
})

test_that("unnamed layers fall back to numbered labels", {
  fit <- fit_of("parallel", 2L)
  expect_equal(names(plot_cluster_omic_profile(fit)), c("Layer 1", "Layer 2"))
})

# ---- the returned plots ----------------------------------------------------

test_that("both renderings return one ggplot per layer for every model type", {
  for (spec in list(list("early", 1L), list("parallel", 2L), list("serial", 2L))) {
    fit <- fit_of(spec[[1]], spec[[2]])
    for (type in c("heatmap", "bar")) {
      p <- plot_cluster_omic_profile(fit, type = type)
      expect_type(p, "list")
      expect_length(p, if (spec[[1]] == "early") 1L else 2L)
      for (g in p) expect_s3_class(g, "ggplot")
    }
  }
})

test_that("each plot carries the data it drew", {
  fit <- fit_of("early", 1L, M = 8)
  p <- plot_cluster_omic_profile(fit, top_n = 5)
  df <- attr(p[[1]], "profile_data")

  expect_s3_class(df, "data.frame")
  expect_named(df, c("feature", "cluster", "value", "mean", "sd", "score"))
  expect_equal(nlevels(df$feature), 5L)          # top_n features
  expect_equal(nrow(df), 5L * nlevels(df$cluster))
  expect_false(anyNA(df$value))
})

test_that("top_n larger than the feature count shows all features", {
  fit <- fit_of("early", 1L, M = 4)
  p <- plot_cluster_omic_profile(fit, top_n = 50)
  expect_equal(nlevels(attr(p[[1]], "profile_data")$feature), 4L)
})

test_that("features are ordered by the importance score", {
  fit <- fit_of("early", 1L, M = 8)
  df <- attr(plot_cluster_omic_profile(fit)[[1]], "profile_data")
  by_feature <- unique(df[, c("feature", "score")])
  # y-axis levels are reversed so the most important sits at the top
  ordered <- by_feature[match(rev(levels(by_feature$feature)), by_feature$feature), ]
  expect_false(is.unsorted(rev(ordered$score)))
})

test_that("the heatmap fill is a per-feature z-score when scale = TRUE", {
  fit <- fit_of("early", 1L, M = 6)
  df <- attr(plot_cluster_omic_profile(fit, scale = TRUE)[[1]], "profile_data")
  by_feat <- split(df$value, df$feature)
  for (v in by_feat) expect_equal(mean(v), 0, tolerance = 1e-8)

  raw <- attr(plot_cluster_omic_profile(fit, scale = FALSE)[[1]], "profile_data")
  expect_equal(raw$value, raw$mean)
})

test_that("the K = 2 fill is centred rather than z-scored by default", {
  # A z-score over two values is always exactly +/-1/sqrt(2), so at K = 2 every
  # tile would saturate and the fill would carry sign but no magnitude. The
  # default therefore centres without rescaling at K = 2 only.
  fit2 <- fit_of("early", 1L, K = 2, M = 8, seed = 25)
  d2 <- attr(plot_cluster_omic_profile(fit2)[[1]], "profile_data")
  expect_true(isTRUE(attr(d2, "centred")))
  # magnitude varies across features, which is the whole point
  expect_gt(stats::sd(abs(d2$value)), 0)
  # and each feature is still centred, so the palette stays diverging
  for (v in split(d2$value, d2$feature)) expect_equal(mean(v), 0, tolerance = 1e-8)

  # an explicit scale = TRUE is still honoured, degenerate or not
  d2e <- attr(plot_cluster_omic_profile(fit2, scale = TRUE)[[1]], "profile_data")
  expect_false(isTRUE(attr(d2e, "centred")))
  expect_equal(unique(round(abs(d2e$value), 8)), round(1 / sqrt(2), 8))

  # at K = 3 the default is the z-score, as documented
  fit3 <- fit_of("early", 1L, K = 3, M = 8, seed = 26)
  d3 <- attr(plot_cluster_omic_profile(fit3)[[1]], "profile_data")
  expect_false(isTRUE(attr(d3, "centred")))
  for (v in split(d3$value, d3$feature)) expect_equal(stats::sd(v), 1, tolerance = 1e-6)
})

test_that("custom labels and colours are honoured", {
  fit <- fit_of("parallel", 2L)
  p <- plot_cluster_omic_profile(fit, layer_names = c("Alpha", "Beta"),
                                 layer_colors = c("#111111", "#222222"))
  expect_equal(names(p), c("Alpha", "Beta"))
  expect_equal(p[[1]]$labels$subtitle, "Alpha")

  df <- attr(plot_cluster_omic_profile(fit, cluster_labels = c("low", "high"))[[1]],
             "profile_data")
  expect_equal(levels(df$cluster), c("low", "high"))
})

test_that("many clusters render without recycling colours", {
  fit <- fit_of("early", 1L, K = 6, M = 8, seed = 24)
  p <- plot_cluster_omic_profile(fit, type = "bar")
  built <- ggplot2::ggplot_build(p[[1]])
  fills <- unique(built$data[[1]]$fill)
  expect_equal(length(fills), 6L)      # one distinct shade per cluster
})

test_that("invalid input is rejected with a clear message", {
  fit <- fit_of("early", 1L)
  expect_error(plot_cluster_omic_profile(fit, cluster_labels = c("only one")),
               "one entry per cluster")
  expect_error(plot_cluster_omic_profile(fit, layer_names = c("a", "b")),
               "one entry per omics layer")
  expect_error(plot_cluster_omic_profile(fit, top_n = 0), "positive")
  expect_error(plot_cluster_omic_profile(structure(list(), class = "lm")),
               "must be a fitted LUCID model")
})

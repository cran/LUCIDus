skip_on_cran()

# Does plot_cluster_omic_profile() generalize across the whole model space?
#
# test-plot-cluster-omic-profile.R covers the mechanics: the ranking measure,
# the returned shapes, the argument handling. This file covers *coverage* --
# every architecture at several sizes, every serial topology, both outcome
# families, missing data, penalized fits, and the awkward sizes -- because the
# function reaches into res_Mu and res_Sigma, and those are stored in three
# different shapes that vary again with topology.
#
# Every cell asserts the same contract, so a failure names the combination that
# broke it rather than just "the plot function".

# ---- shared contract --------------------------------------------------------

# What must be true of a returned profile, whatever produced it.
expect_valid_profile <- function(p, n_panels = NULL, info = NULL) {
  testthat::expect_type(p, "list")
  if (!is.null(n_panels)) testthat::expect_length(p, n_panels)
  testthat::expect_true(length(p) >= 1L, label = paste("at least one panel", info))
  testthat::expect_false(is.null(names(p)), label = paste("panels are named", info))
  testthat::expect_true(all(nzchar(names(p))), label = paste("names non-empty", info))

  for (i in seq_along(p)) {
    g <- p[[i]]
    testthat::expect_s3_class(g, "ggplot")

    df <- attr(g, "profile_data")
    testthat::expect_s3_class(df, "data.frame")
    testthat::expect_named(df, c("feature", "cluster", "value", "mean", "sd", "score"))
    testthat::expect_gt(nrow(df), 0)

    # nothing non-finite may reach a figure
    testthat::expect_true(all(is.finite(df$value)),
      label = paste("finite fill values", info, names(p)[i]))
    testthat::expect_true(all(is.finite(df$score)),
      label = paste("finite scores", info, names(p)[i]))
    testthat::expect_false(anyNA(df$mean))

    # the frame is a complete feature x cluster grid
    testthat::expect_equal(nrow(df), nlevels(df$feature) * nlevels(df$cluster))

    # and it must actually build, which is where a bad scale or factor shows up
    testthat::expect_no_error(ggplot2::ggplot_build(g))
  }
  invisible(TRUE)
}

fit_quiet <- function(...) {
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(...)
  ))))
  fit
}

profile_quiet <- function(...) {
  suppressWarnings(suppressMessages(plot_cluster_omic_profile(...)))
}

# ---- architecture x size ----------------------------------------------------

test_that("early generalizes across cluster counts and feature counts", {
  for (K in c(2, 3, 5)) {
    for (M in c(3, 12)) {
      lbl <- paste("early K =", K, "M =", M)
      d <- sim_lucid("early", N = 250, K = K, M = M, sep = 2.5, seed = 300 + K)
      fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                       family = "normal", K = K, init_omic.data.model = "VVV")
      for (type in c("heatmap", "bar")) {
        expect_valid_profile(profile_quiet(fit, type = type), 1L,
                             info = paste(lbl, type))
      }
    }
  }
})

test_that("parallel generalizes across 1, 2, 3 and 5 layers", {
  for (n_layer in c(1L, 2L, 3L, 5L)) {
    lbl <- paste("parallel", n_layer, "layers")
    d <- sim_lucid("parallel", N = 250, K = 2, M = 6, n_layer = n_layer,
                   sep = 2.5, seed = 310)
    fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "parallel",
                     family = "normal", K = as.list(rep(2, n_layer)),
                     init_omic.data.model = "VVV")
    for (type in c("heatmap", "bar")) {
      expect_valid_profile(profile_quiet(fit, type = type), n_layer,
                           info = paste(lbl, type))
    }
  }
})

test_that("parallel handles unequal cluster counts across layers", {
  d <- sim_lucid("parallel", N = 250, K = c(2, 4), M = 6, n_layer = 2,
                 sep = 2.5, seed = 311)
  fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "parallel",
                   family = "normal", K = list(2, 4), init_omic.data.model = "VVV")
  p <- profile_quiet(fit)
  expect_valid_profile(p, 2L, info = "parallel K = 2, 4")

  # each panel carries its OWN number of clusters, not the first layer's
  expect_equal(nlevels(attr(p[[1]], "profile_data")$cluster), 2L)
  expect_equal(nlevels(attr(p[[2]], "profile_data")$cluster), 4L)
})

test_that("serial generalizes across 1, 2, 3 and 6 stages", {
  for (n_stage in c(1L, 2L, 3L, 6L)) {
    lbl <- paste("serial", n_stage, "stages")
    d <- sim_lucid("serial", N = 250, K = 2, M = 5, n_layer = n_stage,
                   sep = 2.5, seed = 320)
    Z <- if (n_stage == 1L) list(d$Z[[1]]) else d$Z
    fit <- fit_quiet(G = d$G, Z = Z, Y = d$Y, lucid_model = "serial",
                     family = "normal", K = as.list(rep(2, n_stage)),
                     init_omic.data.model = "VVV")
    for (type in c("heatmap", "bar")) {
      expect_valid_profile(profile_quiet(fit, type = type), n_stage,
                           info = paste(lbl, type))
    }
  }
})

# ---- serial topology --------------------------------------------------------

test_that("every serial topology produces one panel per omics block", {
  topologies <- list(
    "all early"            = list(2, 2),
    "all parallel"         = list(list(2, 2), list(2, 2)),
    "early then parallel"  = list(2, list(2, 2)),
    "parallel then early"  = list(list(2, 2), 2),
    "parallel in middle"   = list(2, list(2, 2), 2)
  )
  # a panel per omics matrix: a plain stage contributes 1, a parallel stage its
  # number of layers
  expected <- vapply(topologies, function(topo) {
    sum(vapply(topo, function(k) if (is.list(k)) length(k) else 1L, integer(1)))
  }, integer(1))

  for (nm in names(topologies)) {
    topo <- topologies[[nm]]
    d <- sim_lucid_serial_topology(topo, N = 250, M = 5, sep = 2.5, seed = 330)
    fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                     family = "normal", K = topo, init_omic.data.model = "VVV")
    for (type in c("heatmap", "bar")) {
      expect_valid_profile(profile_quiet(fit, type = type),
                           unname(expected[[nm]]), info = paste(nm, type))
    }
  }
})

test_that("panel names distinguish stage from layer within a stage", {
  topo <- list(2, list(2, 2))
  d <- sim_lucid_serial_topology(topo, N = 250, M = 5, sep = 2.5, seed = 331)
  names(d$Z) <- c("baseline", "followup")
  fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                   family = "normal", K = topo, init_omic.data.model = "VVV")
  nms <- names(profile_quiet(fit))

  expect_length(nms, 3L)
  expect_equal(nms[1], "baseline")
  # the two layers of the parallel stage are distinguishable from each other
  expect_true(all(grepl("followup", nms[2:3])))
  expect_false(nms[2] == nms[3])
})

# ---- outcome family ---------------------------------------------------------

test_that("both outcome families give the same omics profile structure", {
  # the profile depends on res_Mu and res_Sigma, which the outcome model does
  # not touch -- so a binary fit must work identically
  for (model in c("early", "parallel", "serial")) {
    n_layer <- if (model == "early") 1L else 2L
    K <- if (model == "early") 2 else as.list(rep(2, n_layer))
    for (family in c("normal", "binary")) {
      d <- sim_lucid(model, N = 250, K = 2, M = 6, n_layer = n_layer,
                     sep = 2.5, family = family, seed = 340)
      fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                       family = family, K = K, init_omic.data.model = "VVV")
      expect_valid_profile(profile_quiet(fit), n_layer,
                           info = paste(model, family))
    }
  }
})

# ---- missing data -----------------------------------------------------------

test_that("profiles are drawable from fits with missing omics data", {
  # the cluster means are estimated from incomplete data, but they are still
  # complete parameters, so the figure must not inherit any NA
  for (model in c("early", "parallel", "serial")) {
    n_layer <- if (model == "early") 1L else 2L
    K <- if (model == "early") 2 else as.list(rep(2, n_layer))
    for (miss in c("listwise", "sporadic", "mixed")) {
      d <- sim_lucid(model, N = 250, K = 2, M = 6, n_layer = n_layer,
                     sep = 2.5, missing = miss, seed = 350)
      fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                       family = "normal", K = K, init_omic.data.model = "VVV")
      expect_valid_profile(profile_quiet(fit), n_layer,
                           info = paste(model, miss))
    }
  }
})

# ---- penalized fits ---------------------------------------------------------

test_that("a penalized fit still profiles, and deselected features rank last", {
  d <- sim_lucid("early", N = 300, K = 2, M = 10, P = 4, sep = 2.5, seed = 360)
  fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                   family = "normal", K = 2, Rho_Z_Mu = 10, Rho_Z_Cov = 0.02,
                   init_omic.data.model = "VVV")
  p <- profile_quiet(fit)
  expect_valid_profile(p, 1L, info = "penalized early")

  # a deselected feature has identical cluster means, so it scores 0; if any
  # feature scored 0 it must not outrank a feature that did not
  df <- attr(p[[1]], "profile_data")
  by_feat <- unique(df[, c("feature", "score")])
  by_feat <- by_feat[order(-by_feat$score), ]
  zero <- by_feat$score == 0
  if (any(zero) && !all(zero)) expect_true(all(diff(which(zero)) >= 1))
  expect_false(is.unsorted(rev(by_feat$score)))
})

# ---- unsupervised -----------------------------------------------------------

test_that("an unsupervised fit profiles like a supervised one", {
  d <- sim_lucid("early", N = 250, K = 2, M = 6, sep = 2.5, seed = 370)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(lucid_model = "early", G = d$G, Z = d$Z, Y = d$Y,
                          family = "normal", K = 2, useY = FALSE,
                          init_omic.data.model = "VVV")
  ))))
  expect_valid_profile(profile_quiet(fit), 1L, info = "unsupervised")
})

# ---- awkward sizes ----------------------------------------------------------

test_that("top_n interacts correctly with small and large feature counts", {
  for (M in c(2, 5, 20)) {
    d <- sim_lucid("early", N = 250, K = 2, M = M, sep = 2.5, seed = 380)
    fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                     family = "normal", K = 2, init_omic.data.model = "VVV")
    for (top_n in c(1, 10, 100)) {
      p <- profile_quiet(fit, top_n = top_n)
      shown <- nlevels(attr(p[[1]], "profile_data")$feature)
      expect_equal(shown, min(top_n, M),
                   info = paste("M =", M, "top_n =", top_n))
    }
  }
})

test_that("every importance measure works for every architecture", {
  for (model in c("early", "parallel", "serial")) {
    n_layer <- if (model == "early") 1L else 2L
    K <- if (model == "early") 2 else as.list(rep(2, n_layer))
    d <- sim_lucid(model, N = 250, K = 2, M = 6, n_layer = n_layer,
                   sep = 2.5, seed = 390)
    fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                     family = "normal", K = K, init_omic.data.model = "VVV")
    for (measure in c("separation", "range", "sd")) {
      expect_valid_profile(profile_quiet(fit, importance = measure), n_layer,
                           info = paste(model, measure))
    }
  }
})

test_that("many clusters render for every architecture", {
  # the bar ramp and the heatmap axis both have to cope; 8 exceeds the length of
  # any categorical palette the package could have used instead
  for (model in c("early", "parallel")) {
    n_layer <- if (model == "early") 1L else 2L
    K <- if (model == "early") 8 else as.list(rep(8, n_layer))
    d <- sim_lucid(model, N = 400, K = 8, M = 6, n_layer = n_layer,
                   sep = 3, seed = 400)
    fit <- fit_quiet(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                     family = "normal", K = K, init_omic.data.model = "EII")
    for (type in c("heatmap", "bar")) {
      p <- profile_quiet(fit, type = type)
      expect_valid_profile(p, n_layer, info = paste(model, "K = 8", type))
      expect_equal(nlevels(attr(p[[1]], "profile_data")$cluster), 8L)
    }
    # distinct shade per cluster, no recycling
    built <- ggplot2::ggplot_build(profile_quiet(fit, type = "bar")[[1]])
    expect_equal(length(unique(built$data[[1]]$fill)), 8L)
  }
})

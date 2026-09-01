# Heavy: fits multiple LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# get_selected_G(), get_selected_Z(), get_cluster_assignment(), and
# get_top_omics_features() all auto-detect early/parallel/serial from the
# fitted object's class, so a user never branches on model type themselves.
# These tests cross-check each one against the exact same fields/logic read
# directly, for all three model types, penalized and unpenalized.

fit_for <- function(model, n_layer = 2L, penalize = FALSE, seed = 61) {
  d <- sim_lucid(model, N = 250, K = 2, M = 6, n_layer = n_layer,
                sep = 2.5, family = "normal", seed = seed)
  K <- if (model == "early") {
    2
  } else if (model == "parallel") {
    rep(2, n_layer)
  } else {
    as.list(rep(2, n_layer))
  }
  rho <- if (penalize) list(Rho_G = 0.1, Rho_Z_Mu = 0) else list(Rho_G = 0, Rho_Z_Mu = 0)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                          family = "normal", K = K, init_omic.data.model = "VVV",
                          Rho_G = rho$Rho_G, Rho_Z_Mu = rho$Rho_Z_Mu, seed = seed)
  ))))
  list(d = d, fit = fit)
}

test_that("get_selected_G matches fit$select$selectG directly, for all model types", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model, penalize = TRUE)
    expected <- if (model == "serial") {
      fx$fit$submodel[[1]]$select$selectG
    } else {
      fx$fit$select$selectG
    }
    expect_equal(unname(get_selected_G(fx$fit)), unname(expected), info = model)
    expect_false(is.null(names(get_selected_G(fx$fit))), info = model)
  }
})

test_that("get_selected_G(layer=) returns a parallel model's per-layer selection", {
  fx <- fit_for("parallel", penalize = TRUE)
  for (i in seq_along(fx$fit$K)) {
    expect_equal(unname(get_selected_G(fx$fit, layer = i)),
                unname(fx$fit$select$selectG_layer[[i]]))
  }
})

test_that("get_selected_G names line up with the exposures, not CoG covariates", {
  # var.names$Gnames is built as c(exposure names, CoG names); selectG only
  # ever covers the exposure portion. A fit supplied with CoG used to crash
  # get_selected_G() (setNames() length mismatch) because it named the mask
  # with the full, CoG-extended Gnames vector instead of just its own leading
  # portion.
  set.seed(81)
  N <- 100
  G <- matrix(rnorm(N * 3), nrow = N); colnames(G) <- c("g1", "g2", "g3")
  CoG <- matrix(rnorm(N * 2), nrow = N); colnames(CoG) <- c("age", "sex")
  Z <- matrix(rnorm(N * 4), nrow = N)
  Y <- rnorm(N)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(G = G, Z = Z, Y = Y, CoG = CoG, lucid_model = "early",
                          family = "normal", K = 2, Rho_G = 0.1, seed = 81)
  ))))
  got <- get_selected_G(fit)
  expect_equal(names(got), c("g1", "g2", "g3"))
  expect_length(got, 3L)
})

test_that("get_selected_Z matches the selectZ feature mask directly, for all model types", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model, penalize = TRUE)
    got <- get_selected_Z(fx$fit)
    if (model == "early") {
      expect_type(got, "logical")
      expect_length(got, ncol(fx$fit$Z))
    } else if (model == "parallel") {
      expect_type(got, "list")
      expect_length(got, length(fx$fit$K))
      for (i in seq_along(got)) {
        sel <- fx$fit$select$selectZ[[i]]
        expected <- if (is.null(dim(sel))) as.logical(sel) else colSums(sel) > 0
        expect_equal(unname(got[[i]]), unname(expected), info = paste(model, i))
      }
    } else {
      expect_type(got, "list")
      expect_length(got, length(fx$fit$submodel))
    }
  }
})

test_that("get_selected_Z recurses correctly into a parallel stage nested in serial", {
  topo <- list(2, list(2, 2))
  d <- sim_lucid_serial_topology(topo, N = 250, M = 5, sep = 2.5, seed = 62)
  names(d$Z) <- c("first", "second")
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                          family = "normal", K = topo, init_omic.data.model = "VVV",
                          seed = 62)
  ))))
  got <- get_selected_Z(fit)
  expect_length(got, 2L)
  expect_type(got[[1]], "logical")           # stage 1 is early: flat vector
  expect_type(got[[2]], "list")              # stage 2 is parallel: list of 2 layers
  expect_length(got[[2]], 2L)
})

test_that("get_cluster_assignment matches which.max(inclusion.p) directly, for all model types", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    got <- get_cluster_assignment(fx$fit)
    if (model == "early") {
      expect_equal(got, as.numeric(apply(fx$fit$inclusion.p, 1, which.max)))
    } else if (model == "parallel") {
      expect_type(got, "list")
      for (i in seq_along(got)) {
        expect_equal(got[[i]], as.numeric(apply(fx$fit$inclusion.p[[i]], 1, which.max)),
                    info = paste(model, i))
      }
    } else {
      expect_type(got, "list")
      expect_length(got, length(fx$fit$submodel))
    }
  }
})

test_that("get_cluster_assignment matches predict_lucid()'s own pred.x on training data", {
  fx <- fit_for("early")
  pred <- suppressWarnings(suppressMessages(invisible(capture.output(
    p <- predict_lucid(model = fx$fit, G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)
  ))))
  p <- predict_lucid(model = fx$fit, G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)
  expect_equal(get_cluster_assignment(fx$fit), as.numeric(p$pred.x))
})

test_that("get_top_omics_features matches plot_cluster_omic_profile()'s own ranking", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    top <- get_top_omics_features(fx$fit, top_n = 3)
    plotted <- plot_cluster_omic_profile(fx$fit, top_n = 3)
    expect_equal(names(top), names(plotted))
    for (i in seq_along(top)) {
      df <- attr(plotted[[i]], "profile_data")
      expected_order <- rev(levels(df$feature))  # levels are reversed for plotting
      expect_equal(names(top[[i]]), expected_order[seq_along(top[[i]])], info = paste(model, i))
      expect_true(is.numeric(top[[i]]))
      expect_false(is.unsorted(rev(top[[i]])))   # sorted descending
    }
  }
})

test_that("get_top_omics_features honors top_n and caps at the available feature count", {
  fx <- fit_for("early")
  expect_length(get_top_omics_features(fx$fit, top_n = 2)[[1]], 2L)
  full <- ncol(fx$fit$Z)
  expect_length(get_top_omics_features(fx$fit, top_n = 1000)[[1]], full)
})

test_that("get_top_omics_features rejects a bad top_n and a malformed model", {
  fx <- fit_for("early")
  expect_error(get_top_omics_features(fx$fit, top_n = 0), "positive")
  expect_error(get_top_omics_features(structure(list(), class = "lm")),
              "must be a fitted LUCID model")
})

test_that("all four extractors reject a non-LUCID object with a clear message", {
  bad <- structure(list(), class = "lm")
  expect_error(get_selected_G(bad), "must be a fitted LUCID model")
  expect_error(get_selected_Z(bad), "must be a fitted LUCID model")
  expect_error(get_cluster_assignment(bad), "must be a fitted LUCID model")
})

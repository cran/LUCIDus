skip_on_cran()

# predict_lucid() has four independent switches -- g_computation, response,
# whether Y is supplied, and whether Z is supplied -- and three architectures to
# cross them with. The combinations are not interchangeable: g-computation drops
# the Z and Y terms from the E-step entirely, and Z = NULL is legal only in that
# mode.

N_PRED <- 300L
M_PRED <- 4L

pred_fixture <- function(model, n_layer = 2L, family = "normal", seed = 301) {
  d <- sim_lucid(model, N = N_PRED, K = 2, M = M_PRED, n_layer = n_layer,
                 family = family, seed = seed)
  K <- if (model == "early") 2 else as.list(rep(2, n_layer))
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                 family = family, K = K, init_omic.data.model = NULL)
  ))))
  list(d = d, fit = fit)
}

# every predicted cluster label must be a valid 1..K index
expect_labels_in_range <- function(pred_x, K, info = NULL) {
  xs <- if (is.list(pred_x)) pred_x else list(pred_x)
  for (x in xs) {
    x <- as.numeric(unlist(x))
    testthat::expect_true(all(x >= 1 & x <= K),
      label = paste0("cluster labels within 1..", K,
                     if (!is.null(info)) paste0(" [", info, "]")))
    testthat::expect_true(all(x == round(x)), label = "cluster labels are integers")
  }
}

test_that("supervised and unsupervised prediction works for every architecture", {
  for (model in c("early", "parallel", "serial")) {
    for (family in c("normal", "binary")) {
      fx <- pred_fixture(model, family = family)
      lbl <- paste(model, family)

      # Y supplied -> the outcome informs the posterior; Y omitted -> it does not.
      with_y <- predict_lucid(model = fx$fit, lucid_model = model,
                              G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)
      without_y <- predict_lucid(model = fx$fit, lucid_model = model,
                                 G = fx$d$G, Z = fx$d$Z)

      for (p in list(with_y, without_y)) {
        expect_true(all(c("inclusion.p", "pred.x", "pred.y") %in% names(p)), info = lbl)
        expect_labels_in_range(p$pred.x, 2L, info = lbl)
        expect_false(anyNA(unlist(p$pred.y)))
      }
      # supplying Y must actually change something, or useY is being ignored
      expect_false(isTRUE(all.equal(with_y$inclusion.p, without_y$inclusion.p)), info = lbl)
    }
  }
})

test_that("binary response switches between labels and probabilities", {
  for (model in c("early", "parallel", "serial")) {
    fx <- pred_fixture(model, family = "binary")
    labels <- predict_lucid(model = fx$fit, lucid_model = model,
                            G = fx$d$G, Z = fx$d$Z, response = TRUE)
    probs <- predict_lucid(model = fx$fit, lucid_model = model,
                           G = fx$d$G, Z = fx$d$Z, response = FALSE)

    expect_true(all(unlist(labels$pred.y) %in% c(0, 1)), info = model)
    pv <- as.numeric(unlist(probs$pred.y))
    expect_true(all(pv >= 0 & pv <= 1), info = model)
    # probabilities must not already be degenerate 0/1, or the switch is a no-op
    expect_true(any(pv > 0 & pv < 1), info = model)
  }
})

test_that("response has no effect on a normal outcome", {
  # The package gates on is_binary_family() && response, so for a continuous
  # outcome the flag is inert. Asserted so a future change to that gate cannot
  # silently start rounding continuous predictions.
  fx <- pred_fixture("early", family = "normal")
  a <- predict_lucid(model = fx$fit, lucid_model = "early",
                     G = fx$d$G, Z = fx$d$Z, response = TRUE)
  b <- predict_lucid(model = fx$fit, lucid_model = "early",
                     G = fx$d$G, Z = fx$d$Z, response = FALSE)
  expect_equal(a$pred.y, b$pred.y)
})

test_that("g-computation predicts from exposures alone and returns pred.z", {
  for (model in c("early", "parallel", "serial")) {
    fx <- pred_fixture(model)
    g <- predict_lucid(model = fx$fit, lucid_model = model,
                       G = fx$d$G, Z = NULL, g_computation = TRUE)

    expect_true("pred.z" %in% names(g), info = model)
    expect_false(is.null(g$pred.z), info = model)
    expect_labels_in_range(g$pred.x, 2L, info = model)
    expect_false(anyNA(unlist(g$pred.z)), info = model)
    expect_false(anyNA(unlist(g$inclusion.p)), info = model)

    # Z and Y are ignored in this mode, so passing them changes nothing.
    g2 <- predict_lucid(model = fx$fit, lucid_model = model, G = fx$d$G,
                        Z = fx$d$Z, Y = fx$d$Y, g_computation = TRUE)
    expect_equal(g$inclusion.p, g2$inclusion.p, info = model)
  }
})

test_that("counterfactual exposures move the g-computation prediction", {
  # The point of g-computation is that changing G changes the predicted omics
  # and outcome. A test that only checked it runs would pass on a constant.
  fx <- pred_fixture("early")
  hi <- fx$d$G; hi[, 1] <- hi[, 1] + 2
  lo <- fx$d$G; lo[, 1] <- lo[, 1] - 2

  p_hi <- predict_lucid(model = fx$fit, lucid_model = "early", G = hi,
                        Z = NULL, g_computation = TRUE)
  p_lo <- predict_lucid(model = fx$fit, lucid_model = "early", G = lo,
                        Z = NULL, g_computation = TRUE)

  expect_false(isTRUE(all.equal(p_hi$inclusion.p, p_lo$inclusion.p)))
  expect_false(isTRUE(all.equal(p_hi$pred.y, p_lo$pred.y)))
  # the first exposure carries the generating signal, so raising it must shift
  # cluster membership in a consistent direction
  expect_true(mean(p_hi$inclusion.p[, 2]) != mean(p_lo$inclusion.p[, 2]))
})

test_that("non-g-computation prediction still requires Z", {
  for (model in c("early", "parallel", "serial")) {
    fx <- pred_fixture(model)
    expect_error(
      predict_lucid(model = fx$fit, lucid_model = model, G = fx$d$G, Z = NULL),
      info = model
    )
  }
})

test_that("predicting on the training data reproduces the fitted assignment", {
  # Extracting cluster membership by predicting on the training data is a
  # documented use, so it must agree with what the fit already reported.
  for (model in c("early", "parallel")) {
    fx <- pred_fixture(model)
    p <- predict_lucid(model = fx$fit, lucid_model = model,
                       G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)
    if (model == "early") {
      expect_equal(as.numeric(p$pred.x), as.numeric(max.col(fx$fit$inclusion.p)))
    } else {
      for (i in seq_along(p$pred.x)) {
        expect_equal(as.numeric(p$pred.x[[i]]),
                     as.numeric(max.col(as.matrix(fx$fit$inclusion.p[[i]]))))
      }
    }
  }
})

test_that("prediction tolerates missing omics in the new data", {
  for (model in c("early", "parallel")) {
    n_layer <- if (model == "early") 1L else 2L
    d <- sim_lucid(model, N = N_PRED, K = 2, M = M_PRED, n_layer = n_layer, seed = 302)
    K <- if (model == "early") 2 else as.list(rep(2, n_layer))
    suppressWarnings(suppressMessages(invisible(capture.output(
      fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                   family = "normal", K = K, init_omic.data.model = NULL)))))

    dm <- sim_lucid(model, N = N_PRED, K = 2, M = M_PRED, n_layer = n_layer,
                    missing = "mixed", seed = 302)
    p <- predict_lucid(model = fit, lucid_model = model, G = dm$G, Z = dm$Z)

    expect_false(anyNA(unlist(p$inclusion.p)), info = model)
    expect_labels_in_range(p$pred.x, 2L, info = model)
  }
})

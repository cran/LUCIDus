# Heavy: fits multiple LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# What predict_lucid() accepts, and what it says when it declines.
#
# The rule is one sentence: Z is required for every model type, and only
# g_computation = TRUE relaxes it. Y is always optional -- omitting it switches
# the posterior from supervised to unsupervised.
#
# The enforcement used to differ by model type. Early and parallel reported a
# missing Z correctly; serial compared length(Z) against length(K) first, and
# length(NULL) is 0, so a forgotten argument was reported as a topology
# mismatch and sent the user to check their K. These tests pin the rule and the
# uniformity of the message.

N_ARG <- 250L

fit_for <- function(model, n_layer = 2L, seed = 7) {
  d <- sim_lucid(model, N = N_ARG, K = 2, M = 4, n_layer = n_layer,
                 sep = 2.5, family = "normal", seed = seed)
  K <- if (model == "early") 2 else as.list(rep(2, n_layer))
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                 family = "normal", K = K, init_omic.data.model = "VVV")
  ))))
  list(d = d, fit = fit)
}

quiet_predict <- function(...) {
  suppressWarnings(suppressMessages(invisible(capture.output(
    out <- predict_lucid(...)
  ))))
  out
}

test_that("Z plus Y works for every model type", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    p <- quiet_predict(model = fx$fit, lucid_model = model,
                       G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)
    expect_true(all(c("inclusion.p", "pred.x", "pred.y") %in% names(p)), info = model)
    expect_null(p$pred.z, info = model)
  }
})

test_that("omitting Y is always legal and gives unsupervised prediction", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    with_y <- quiet_predict(model = fx$fit, lucid_model = model,
                            G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)
    no_y <- quiet_predict(model = fx$fit, lucid_model = model,
                          G = fx$d$G, Z = fx$d$Z)

    # pred.y is still produced without Y: the model predicts the outcome from
    # the cluster assignment it just made.
    expect_false(is.null(no_y$pred.y), info = model)
    expect_false(anyNA(unlist(no_y$pred.y)), info = model)

    # and supplying Y must actually change the posterior, or useY is ignored
    expect_false(isTRUE(all.equal(with_y$inclusion.p, no_y$inclusion.p)), info = model)
  }
})

test_that("omitting Z errors with the same message for every model type", {
  msgs <- character(0)
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    e <- tryCatch(
      quiet_predict(model = fx$fit, lucid_model = model, G = fx$d$G, Z = NULL, Y = fx$d$Y),
      error = function(e) conditionMessage(e)
    )
    expect_true(is.character(e), info = model)
    expect_match(e, "'Z' is required", info = model)
    # the message must name the one mode that relaxes the rule, so the error
    # teaches the rule rather than only reporting the failure
    expect_match(e, "g_computation", info = model)
    msgs <- c(msgs, e)
  }
  # serial used to say "Z and K should be two lists of the same length"
  expect_length(unique(msgs), 1L)
})

test_that("a real Z/K length mismatch is still reported as one", {
  # The reordering must not mask a genuine topology error behind the new guard.
  fx <- fit_for("serial", n_layer = 2L)
  too_few <- fx$d$Z[1]                      # 1 element for a 2-stage model
  expect_error(
    quiet_predict(model = fx$fit, lucid_model = "serial",
                  G = fx$d$G, Z = too_few, Y = fx$d$Y),
    "same length"
  )
})

test_that("a non-list Z for serial is reported as a type error, not a length one", {
  fx <- fit_for("serial", n_layer = 2L)
  expect_error(
    quiet_predict(model = fx$fit, lucid_model = "serial",
                  G = fx$d$G, Z = fx$d$Z[[1]], Y = fx$d$Y),
    "should be a list"
  )
})

test_that("g-computation is the one mode that accepts Z = NULL", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    p <- quiet_predict(model = fx$fit, lucid_model = model,
                       G = fx$d$G, Z = NULL, Y = NULL, g_computation = TRUE)
    expect_true(all(c("inclusion.p", "pred.x", "pred.y", "pred.z") %in% names(p)),
                info = model)
    expect_false(anyNA(unlist(p$pred.z)), info = model)
  }
})

test_that("Z and Y are genuinely ignored under g-computation", {
  # The function prints a notice saying they will not be used. Assert that the
  # notice is telling the truth, so nobody reads a difference into results from
  # a call that passed them.
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    bare <- quiet_predict(model = fx$fit, lucid_model = model,
                          G = fx$d$G, Z = NULL, Y = NULL, g_computation = TRUE)
    loaded <- quiet_predict(model = fx$fit, lucid_model = model,
                            G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y, g_computation = TRUE)
    expect_equal(bare$inclusion.p, loaded$inclusion.p, info = model)
    expect_equal(bare$pred.y, loaded$pred.y, info = model)
  }
})

test_that("a single-stage serial model is declined with an actionable message", {
  # It fits without complaint, but the stage loop cannot predict it: the
  # branches are first / middle / last and one stage is both first and last, so
  # pred.y is never assigned. Rather than surface that as
  # "object 'pred.y' not found", the model type is declined and the user is
  # pointed at the equivalent early or parallel fit.
  d <- sim_lucid("serial", N = N_ARG, K = 2, M = 4, n_layer = 1,
                 sep = 2.5, seed = 7)
  Z1 <- list(d$Z[[1]])
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit1 <- lucid(G = d$G, Z = Z1, Y = d$Y, lucid_model = "serial",
                  family = "normal", K = list(2), init_omic.data.model = "VVV")
  ))))
  expect_s3_class(fit1, "lucid_serial")     # it really does fit

  expect_error(
    quiet_predict(model = fit1, lucid_model = "serial", G = d$G, Z = Z1, Y = d$Y),
    "at least two stages"
  )
})

test_that("predict_lucid() rejects a wrong-class model for parallel, matching early", {
  # Early already checked inherits(model, "early_lucid"); parallel had no
  # analogous check and would fail deep inside the prediction machinery on a
  # mismatched object instead of with a clear message naming the problem.
  fx <- fit_for("parallel")
  bad_model <- fx$fit
  class(bad_model) <- "early_lucid"
  expect_error(
    quiet_predict(model = bad_model, lucid_model = "parallel",
                 G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y),
    "lucid_parallel"
  )
})

test_that("predict_lucid() rejects a wrong-width Z for early, matching parallel", {
  # Parallel already checked length(Z) against the model's layer count; early
  # never checked ncol(Z) against the fitted model's omics width, so a
  # wrong-width Z used to fail inside Estep_early() with a linear-algebra
  # error instead of a clear one naming the expected width.
  fx <- fit_for("early", n_layer = 1L)
  too_narrow <- fx$d$Z[, 1, drop = FALSE]
  expect_error(
    quiet_predict(model = fx$fit, lucid_model = "early",
                 G = fx$d$G, Z = too_narrow, Y = fx$d$Y),
    "should have"
  )
})

test_that("the legality rule holds as a whole matrix", {
  # The table the vignettes document, asserted directly.
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    ok <- function(expr) !inherits(try(suppressWarnings(suppressMessages(
      invisible(capture.output(force(expr))))), silent = TRUE), "try-error")

    expect_true(ok(predict_lucid(model = fx$fit, lucid_model = model,
                                 G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)), info = model)
    expect_true(ok(predict_lucid(model = fx$fit, lucid_model = model,
                                 G = fx$d$G, Z = fx$d$Z)), info = model)
    expect_false(ok(predict_lucid(model = fx$fit, lucid_model = model,
                                  G = fx$d$G, Z = NULL, Y = fx$d$Y)), info = model)
    expect_false(ok(predict_lucid(model = fx$fit, lucid_model = model,
                                  G = fx$d$G, Z = NULL)), info = model)
    expect_true(ok(predict_lucid(model = fx$fit, lucid_model = model,
                                 G = fx$d$G, Z = NULL, Y = NULL,
                                 g_computation = TRUE)), info = model)
  }
})

test_that("lucid_model is auto-detected from model's class when omitted", {
  for (model in c("early", "parallel", "serial")) {
    fx <- fit_for(model)
    with_arg <- quiet_predict(model = fx$fit, lucid_model = model,
                              G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)
    without_arg <- quiet_predict(model = fx$fit,
                                 G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y)
    expect_equal(with_arg, without_arg, info = model)
  }
})

test_that("an explicit, wrong lucid_model is still rejected when supplied", {
  fx <- fit_for("early")
  expect_error(
    quiet_predict(model = fx$fit, lucid_model = "parallel",
                 G = fx$d$G, Z = fx$d$Z, Y = fx$d$Y),
    "lucid_parallel"
  )
})

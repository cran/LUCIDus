# Note: predict_lucid() now reports a missing Z with one message for every
# model type -- "Input data 'Z' is required for prediction..." -- naming
# g_computation as the one mode that relaxes the rule. These expectations
# match on the stable part of that message.

test_that("early g-computation accepts Z = NULL and ignores Z/Y inputs", {
  G <- sim_data$G[1:120, ]
  Z <- sim_data$Z[1:120, ]
  Y <- sim_data$Y_normal[1:120, ]
  CoY <- sim_data$Covariate[1:120, ]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G,
      Z = Z,
      Y = Y,
      CoY = CoY,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  invisible(capture.output(
    pred_with_zy <- predict_lucid(
      model = fit,
      lucid_model = "early",
      G = G,
      Z = Z,
      Y = Y,
      CoY = CoY,
      g_computation = TRUE,
      response = FALSE
    )
  ))

  pred_no_zy <- predict_lucid(
    model = fit,
    lucid_model = "early",
    G = G,
    Z = NULL,
    Y = NULL,
    CoY = CoY,
    g_computation = TRUE,
    response = FALSE
  )

  expect_equal(pred_with_zy$inclusion.p, pred_no_zy$inclusion.p, tolerance = 1e-8)
  expect_equal(pred_with_zy$pred.y, pred_no_zy$pred.y, tolerance = 1e-8)
  expect_equal(pred_with_zy$pred.z, pred_no_zy$pred.z, tolerance = 1e-8)
  expect_equal(nrow(pred_no_zy$pred.z), nrow(G))
  expect_equal(ncol(pred_no_zy$pred.z), ncol(fit$res_Mu))
  expect_equal(rowSums(pred_no_zy$inclusion.p), rep(1, nrow(G)), tolerance = 1e-6)
})

test_that("parallel g-computation accepts Z = NULL and returns layer-wise pred.z", {
  set.seed(2026)
  N <- 60
  G <- matrix(rnorm(N * 3), nrow = N)
  Z <- list(
    matrix(rnorm(N * 5), nrow = N),
    matrix(rnorm(N * 4), nrow = N)
  )
  Y <- rnorm(N)
  CoY <- matrix(rnorm(N * 2), nrow = N)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = Z,
      Y = Y,
      CoY = CoY,
      family = "normal",
      K = c(2, 2),
      seed = 2026,
      max_itr = 8,
      tol = 1e-1
    )
  )))

  invisible(capture.output(
    pred_with_zy <- predict_lucid(
      model = fit,
      lucid_model = "parallel",
      G = G,
      Z = Z,
      Y = Y,
      CoY = CoY,
      g_computation = TRUE,
      response = FALSE
    )
  ))

  pred_no_zy <- predict_lucid(
    model = fit,
    lucid_model = "parallel",
    G = G,
    Z = NULL,
    Y = NULL,
    CoY = CoY,
    g_computation = TRUE,
    response = FALSE
  )

  expect_equal(pred_with_zy$pred.y, pred_no_zy$pred.y, tolerance = 1e-8)
  expect_equal(pred_with_zy$inclusion.p, pred_no_zy$inclusion.p, tolerance = 1e-8)
  expect_equal(length(pred_no_zy$pred.z), 2)

  for (i in seq_along(pred_no_zy$pred.z)) {
    expect_equal(pred_with_zy$pred.z[[i]], pred_no_zy$pred.z[[i]], tolerance = 1e-8)
    expect_equal(nrow(pred_no_zy$pred.z[[i]]), N)
    k_i <- fit$K[i]
    expected_z <- if (nrow(fit$res_Mu[[i]]) == k_i) ncol(fit$res_Mu[[i]]) else nrow(fit$res_Mu[[i]])
    expect_equal(ncol(pred_no_zy$pred.z[[i]]), expected_z)
    expect_equal(rowSums(pred_no_zy$inclusion.p[[i]]), rep(1, N), tolerance = 1e-6)
  }

  # Regression: g-comp should also work when Mu is stored as feature x cluster.
  fit_transposed <- fit
  for (i in seq_along(fit_transposed$res_Mu)) {
    mu_t <- t(fit_transposed$res_Mu[[i]])
    rownames(mu_t) <- paste0("feature_", i, "_", seq_len(nrow(mu_t)))
    colnames(mu_t) <- paste0("cluster", seq_len(ncol(mu_t)))
    fit_transposed$res_Mu[[i]] <- mu_t
  }
  pred_transposed <- predict_lucid(
    model = fit_transposed,
    lucid_model = "parallel",
    G = G,
    Z = NULL,
    Y = NULL,
    CoY = CoY,
    g_computation = TRUE,
    response = FALSE
  )
  for (i in seq_along(pred_transposed$pred.z)) {
    expect_equal(ncol(pred_transposed$pred.z[[i]]), nrow(fit_transposed$res_Mu[[i]]))
    expect_equal(colnames(pred_transposed$pred.z[[i]]), rownames(fit_transposed$res_Mu[[i]]))
  }
})

test_that("non-g-computation still requires Z for early and parallel", {
  G <- sim_data$G[1:80, ]
  Z <- sim_data$Z[1:80, ]
  Y <- sim_data$Y_normal[1:80, ]

  suppressWarnings(invisible(capture.output(
    fit_early <- estimate_lucid(
      G = G,
      Z = Z,
      Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  expect_error(
    predict_lucid(
      model = fit_early,
      lucid_model = "early",
      G = G,
      Z = NULL,
      g_computation = FALSE
    ),
    "'Z' is required"
  )

  set.seed(88)
  Gp <- matrix(rnorm(160), nrow = 80)
  Zp <- list(matrix(rnorm(240), nrow = 80), matrix(rnorm(240), nrow = 80))
  Yp <- rnorm(80)

  suppressWarnings(invisible(capture.output(
    fit_parallel <- estimate_lucid(
      lucid_model = "parallel",
      G = Gp,
      Z = Zp,
      Y = Yp,
      family = "normal",
      K = c(2, 2),
      seed = 88,
      max_itr = 6,
      tol = 1e-1
    )
  )))

  expect_error(
    predict_lucid(
      model = fit_parallel,
      lucid_model = "parallel",
      G = Gp,
      Z = NULL,
      g_computation = FALSE
    ),
    "'Z' is required"
  )
})

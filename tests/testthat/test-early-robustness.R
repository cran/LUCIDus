# Additional robustness tests for early LUCID workflow

test_that("early missing data with LOD imputation returns finite posteriors", {
  G <- sim_data$G[1:180, ]
  Z <- sim_data$Z[1:180, ]
  Y <- sim_data$Y_normal[1:180, ]

  # Introduce sporadic missingness
  Z[1:30, 1] <- NA
  Z[31:60, 2] <- NA

  set.seed(1008)
  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      init_impute = "lod",
      useY = TRUE,
      seed = 1008
    )
  )))

  expect_s3_class(fit, "early_lucid")
  expect_true(all(is.finite(fit$inclusion.p)))
  expect_equal(nrow(fit$inclusion.p), nrow(G))
})

test_that("early missing data keeps all-missing rows as NA in stored Z", {
  G <- sim_data$G[1:120, ]
  Z <- sim_data$Z[1:120, ]
  Y <- sim_data$Y_normal[1:120, ]

  # Entire first row missing in omics
  Z[1, ] <- NA

  set.seed(1008)
  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      init_impute = "mix",
      seed = 1008
    )
  )))

  expect_s3_class(fit, "early_lucid")
  expect_true(all(is.na(fit$Z[1, ])))
})

test_that("tune_lucid early over K grid returns one row per K candidate", {
  G <- sim_data$G[1:160, ]
  Z <- sim_data$Z[1:160, ]
  Y <- sim_data$Y_normal[1:160, ]

  suppressWarnings(invisible(capture.output(
    tuned <- tune_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2:3,
      seed = 1008,
      useY = TRUE
    )
  )))

  expect_equal(nrow(tuned$tune_list), 2)
  expect_true("BIC" %in% colnames(tuned$tune_list))
  expect_s3_class(tuned$best_model, "early_lucid")
})

test_that("tune_lucid early over penalty grid returns full combinations", {
  G <- sim_data$G[1:140, ]
  Z <- sim_data$Z[1:140, ]
  Y <- sim_data$Y_normal[1:140, ]

  suppressWarnings(invisible(capture.output(
    tuned <- tune_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      Rho_G = c(0, 0.01),
      Rho_Z_Mu = c(0, 1),
      Rho_Z_Cov = 0,
      seed = 1008,
      useY = TRUE
    )
  )))

  expect_equal(nrow(tuned$tune_list), 4)
  expect_true("BIC" %in% colnames(tuned$tune_list))
  expect_type(tuned$res_model, "list")
})

test_that("lucid wrapper chooses one K from candidate vector", {
  G <- sim_data$G[1:180, ]
  Z <- sim_data$Z[1:180, ]
  Y <- sim_data$Y_normal[1:180, ]

  suppressWarnings(invisible(capture.output(
    fit <- lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2:3,
      seed = 1008
    )
  )))

  expect_s3_class(fit, "early_lucid")
  expect_true(fit$K %in% c(2, 3))
})

test_that("lucid wrapper with Rho_G vector returns logical selectG", {
  G <- sim_data$G[1:180, ]
  Z <- sim_data$Z[1:180, ]
  Y <- sim_data$Y_binary[1:180, ]
  cov <- sim_data$Covariate[1:180, ]

  suppressWarnings(invisible(capture.output(
    fit <- lucid(
      G = G, Z = Z, Y = Y,
      CoY = cov,
      lucid_model = "early",
      family = "binary",
      K = 2,
      Rho_G = c(0, 0.05),
      seed = 1008
    )
  )))

  expect_type(fit$select$selectG, "logical")
  expect_true(length(fit$select$selectG) >= 1)
  expect_true(length(fit$select$selectG) <= ncol(G))
})

test_that("lucid early preserves tuned Rho values in refit model metadata", {
  G <- sim_data$G[1:160, ]
  Z <- sim_data$Z[1:160, ]
  Y <- sim_data$Y_normal[1:160, ]

  suppressWarnings(invisible(capture.output(
    fit <- lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      Rho_G = 0.01,
      Rho_Z_Mu = 0.1,
      Rho_Z_Cov = 0.01,
      max_itr = 8,
      max_tot.itr = 20,
      tol = 1e-2,
      seed = 1008
    )
  )))

  expect_equal(fit$Rho$Rho_G, 0.01)
  expect_equal(fit$Rho$Rho_Z_Mu, 0.1)
  expect_equal(fit$Rho$Rho_Z_Cov, 0.01)
})

test_that("lucid early penalty refit uses tuned scalar K (not candidate grid)", {
  G <- sim_data$G[1:160, ]
  Z <- sim_data$Z[1:160, ]
  Y <- sim_data$Y_normal[1:160, ]

  suppressWarnings(invisible(capture.output(
    fit <- lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2:3,
      Rho_G = c(0, 0.01),
      Rho_Z_Mu = 0,
      Rho_Z_Cov = 0,
      max_itr = 8,
      max_tot.itr = 20,
      tol = 1e-2,
      seed = 1008
    )
  )))

  expect_s3_class(fit, "early_lucid")
  expect_true(is.numeric(fit$K))
  expect_equal(length(fit$K), 1)
  expect_true(fit$K %in% c(2, 3))
})

test_that("lucid wrapper rejects negative penalties", {
  G <- sim_data$G[1:120, ]
  Z <- sim_data$Z[1:120, ]
  Y <- sim_data$Y_normal[1:120, ]

  expect_error(
    lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      Rho_G = -0.01
    ),
    "greater than or equal to 0"
  )
})

test_that("estimate_lucid binary rejects non 0-1 outcomes", {
  G <- sim_data$G[1:120, ]
  Z <- sim_data$Z[1:120, ]
  Y <- sim_data$Y_binary[1:120, ] + 1

  expect_error(
    estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "binary",
      K = 2
    ),
    "coded as 0 and 1|contain only 0s and 1s"
  )
})

test_that("early g-computation returns pred.z and normalized inclusion probabilities", {
  G <- sim_data$G[1:150, ]
  Z <- sim_data$Z[1:150, ]
  Y <- sim_data$Y_binary[1:150, ]
  cov <- sim_data$Covariate[1:150, ]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      CoY = cov,
      lucid_model = "early",
      family = "binary",
      K = 2,
      seed = 1008
    )
  )))

  pred <- predict_lucid(
    model = fit,
    lucid_model = "early",
    G = G,
    Z = Z,
    Y = NULL,
    CoY = cov,
    g_computation = TRUE
  )

  expect_true("pred.z" %in% names(pred))
  expect_equal(nrow(pred$pred.z), nrow(G))
  expect_equal(rowSums(pred$inclusion.p), rep(1, nrow(G)), tolerance = 1e-6)
})

test_that("summary_lucid early exposes top-level and nested BIC consistently", {
  G <- sim_data$G[1:140, ]
  Z <- sim_data$Z[1:140, ]
  Y <- sim_data$Y_normal[1:140, ]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_early")
  expect_true(is.numeric(s$BIC))
  expect_equal(s$BIC, s$model_fit$BIC, tolerance = 1e-8)
})

test_that("summary_lucid early reports listwise and sporadic missing-data profile", {
  G <- sim_data$G[1:120, ]
  Z <- sim_data$Z[1:120, ]
  Y <- sim_data$Y_normal[1:120, ]
  Z[1, ] <- NA
  Z[2:4, 1] <- NA

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_early")
  expect_false(is.null(s$missing_data))
  expect_equal(s$missing_data$listwise_rows, 1)
  expect_equal(s$missing_data$sporadic_rows, 3)
  out <- capture.output(print(s))
  expect_true(any(grepl("Missing-data profile", out)))
})

test_that("estimate_lucid early stores EM control settings for bootstrap reuse", {
  G <- sim_data$G[1:120, ]
  Z <- sim_data$Z[1:120, ]
  Y <- sim_data$Y_normal[1:120, ]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      tol = 1e-2,
      max_itr = 9,
      max_tot.itr = 30,
      seed = 1008
    )
  )))

  expect_true(is.list(fit$em_control))
  expect_equal(fit$em_control$tol, 1e-2)
  expect_equal(fit$em_control$max_itr, 9)
  expect_equal(fit$em_control$max_tot.itr, 30)
})

test_that("summary_lucid early prints intercept (not forced cluster1=0) when CoY is included", {
  G <- sim_data$G[1:140, ]
  Z <- sim_data$Z[1:140, ]
  Y <- sim_data$Y_normal[1:140, ]
  CoY <- sim_data$Covariate[1:140, , drop = FALSE]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      CoY = CoY,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  s <- summary_lucid(fit)
  out <- capture.output(print(s))
  expect_true(any(grepl("^\\(Intercept\\)\\s", out)))
  expect_false(any(grepl("^cluster1\\s+0", out)))
})

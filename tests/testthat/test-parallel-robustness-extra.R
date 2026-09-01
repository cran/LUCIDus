# Heavy: fits multiple LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# Additional robustness tests for parallel LUCID

test_that("check_na for parallel classifies row-level missing patterns correctly", {
  Z1 <- matrix(rnorm(40), nrow = 10)
  Z2 <- matrix(rnorm(40), nrow = 10)

  Z1[1, ] <- NA      # all missing
  Z1[2, 1] <- NA     # partial missing
  Z2[3, ] <- NA      # all missing
  Z2[4, 2] <- NA     # partial missing

  na_pat <- check_na(list(Z1, Z2), lucid_model = "parallel")

  expect_equal(na_pat$indicator_na[[1]][1], 3)
  expect_equal(na_pat$indicator_na[[1]][2], 2)
  expect_equal(na_pat$indicator_na[[1]][5], 1)
  expect_equal(na_pat$indicator_na[[2]][3], 3)
  expect_equal(na_pat$indicator_na[[2]][4], 2)
  expect_true(all(na_pat$impute_flag == c(TRUE, TRUE)))
})

test_that("parallel LOD imputation fills partial rows but leaves all-missing rows NA", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)

  Z1[1, ] <- NA
  Z1[2, 1:2] <- NA

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = c(2, 2),
      family = "normal",
      init_impute = "lod",
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  # D7: a row with no measured omics value at all must stay NA under every
  # initializer.  The LOD path used to fill it with LOD/sqrt(2), returning
  # fabricated data for participants who were never measured.
  expect_true(all(is.na(fit$Z[[1]][1, ])))
  expect_true(all(is.finite(fit$Z[[1]][2, ])))
})

test_that("parallel tune_lucid carries penalty grid and returns model with tuned Rho", {
  set.seed(1008)
  G <- matrix(rnorm(120), nrow = 30)
  Z1 <- matrix(rnorm(180), nrow = 30)
  Z2 <- matrix(rnorm(180), nrow = 30)
  Y <- rnorm(30)

  suppressWarnings(invisible(capture.output(
    tuned <- tune_lucid(
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = list(2:3, 2),
      Rho_G = c(0, 0.05),
      Rho_Z_Mu = c(0, 0.1),
      Rho_Z_Cov = 0,
      max_itr = 6,
      tol = 1e-1,
      seed = 1008,
      useY = TRUE
    )
  )))

  expect_equal(nrow(tuned$tune_K), 8)
  expect_true(all(c("Rho_G", "Rho_Z_Mu", "Rho_Z_Cov", "BIC") %in% colnames(tuned$tune_K)))
  expect_s3_class(tuned$model_opt, "lucid_parallel")
  expect_true(tuned$model_opt$Rho$Rho_G %in% c(0, 0.05))
  expect_true(tuned$model_opt$Rho$Rho_Z_Mu %in% c(0, 0.1))
})

test_that("parallel lucid wrapper tunes K and penalty vectors together", {
  set.seed(1008)
  G <- matrix(rnorm(120), nrow = 30)
  Z1 <- matrix(rnorm(180), nrow = 30)
  Z2 <- matrix(rnorm(180), nrow = 30)
  Y <- rnorm(30)

  suppressWarnings(invisible(capture.output(
    fit <- lucid(
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = list(2:3, 2),
      Rho_G = c(0, 0.05),
      Rho_Z_Mu = 0,
      Rho_Z_Cov = 0,
      max_itr = 6,
      tol = 1e-1,
      seed = 1008,
      useY = TRUE
    )
  )))

  expect_s3_class(fit, "lucid_parallel")
  expect_true(fit$K[1] %in% c(2, 3))
  expect_equal(fit$K[2], 2)
  expect_true(fit$Rho$Rho_G %in% c(0, 0.05))
})

test_that("parallel selection objects are consistent with feature dimensions", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = c(2, 2),
      family = "normal",
      Rho_G = 0.1,
      Rho_Z_Mu = 0.1,
      Rho_Z_Cov = 0,
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  expect_type(fit$select$selectG, "logical")
  expect_true(is.list(fit$select$selectZ))
  expect_true(is.list(fit$select$selectG_layer))
  expect_equal(length(fit$select$selectG_layer), 2)
  expect_equal(length(fit$select$selectZ), 2)
  expect_equal(length(fit$select$selectG), ncol(G))
  expect_equal(length(fit$select$selectG_layer[[1]]), ncol(G))
  expect_equal(length(fit$select$selectG_layer[[2]]), ncol(G))
  expect_equal(dim(fit$select$selectZ[[1]]), c(fit$K[1], ncol(Z1)))
  expect_equal(dim(fit$select$selectZ[[2]]), c(fit$K[2], ncol(Z2)))
})

test_that("parallel summary reports selected-feature tables with valid bounds", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = c(2, 2),
      family = "normal",
      Rho_G = 0.1,
      Rho_Z_Mu = 0.1,
      Rho_Z_Cov = 0,
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_parallel")
  expect_equal(length(s$feature_selection$Z), 2)
  expect_equal(nrow(s$feature_selection$Z[[1]]), ncol(Z1))
  expect_equal(nrow(s$feature_selection$Z[[2]]), ncol(Z2))
  expect_true(all(s$feature_selection$Z[[1]]$Selected_in_clusters >= 0))
  expect_true(all(s$feature_selection$Z[[1]]$Selected_in_clusters <= fit$K[1]))
  expect_true(all(s$feature_selection$Z[[2]]$Selected_in_clusters >= 0))
  expect_true(all(s$feature_selection$Z[[2]]$Selected_in_clusters <= fit$K[2]))
})

test_that("parallel summary reports per-layer listwise and sporadic missing profile", {
  skip_if_not_installed("mix")
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)

  Z1[1, ] <- NA
  Z1[2, 1:2] <- NA
  Z2[3, ] <- NA
  Z2[4:5, 1] <- NA

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = c(2, 2),
      family = "normal",
      init_impute = "mix",
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_parallel")
  expect_false(is.null(s$missing_data))
  expect_true("layer_summary" %in% names(s$missing_data))
  expect_equal(s$missing_data$layer_summary$listwise_rows[1], 1)
  expect_equal(s$missing_data$layer_summary$sporadic_rows[1], 1)
  expect_equal(s$missing_data$layer_summary$listwise_rows[2], 1)
  expect_equal(s$missing_data$layer_summary$sporadic_rows[2], 2)
  out <- capture.output(print(s))
  expect_true(any(grepl("Missing-data profile by layer", out)))
})

# LUCID tuning smoke test for early model (normal)

test_that("tune_lucid and lucid wrapper run for early normal with small grid", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000), nrow = 100)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)
  Y <- rnorm(100)

  suppressWarnings(invisible(capture.output(
    tune_list1 <- tune_lucid(
      G = G, Z = Z1, Y = Y,
      K = 2:3,
      CoG = CoG, CoY = CoY,
      lucid_model = "early",
      family = "normal",
      init_omic.data.model = "VVV",
      seed = i,
      init_impute = "mix",
      init_par = "mclust",
      useY = TRUE
    )
  )))

  suppressWarnings(invisible(capture.output(
    best1 <- lucid(
      G = G, Z = Z1, Y = Y,
      K = 2:3,
      CoG = CoG, CoY = CoY,
      lucid_model = "early",
      family = "normal",
      init_omic.data.model = "VVV",
      seed = i,
      init_impute = "mix",
      init_par = "mclust",
      useY = TRUE
    )
  )))

  expect_s3_class(best1, "early_lucid")
  expect_equal(nrow(tune_list1$tune_list), 2)
  expect_true("BIC" %in% names(tune_list1$tune_list))
  expect_s3_class(tune_list1$best_model, "early_lucid")
})

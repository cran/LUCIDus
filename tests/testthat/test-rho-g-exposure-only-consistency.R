test_that("early Rho_G selection requires >=2 exposures even with multi-column CoG", {
  set.seed(1008)
  N <- 40
  G <- matrix(rnorm(N), nrow = N)            # one exposure only
  CoG <- matrix(rnorm(N * 2), nrow = N)      # extra covariates should not bypass rule
  Z <- matrix(rnorm(N * 4), nrow = N)
  Y <- rnorm(N)

  expect_error(
    estimate_lucid(
      lucid_model = "early",
      G = G,
      Z = Z,
      Y = Y,
      CoG = CoG,
      family = "normal",
      K = 2,
      Rho_G = 0.1,
      max_itr = 6,
      tol = 1e-1,
      seed = 1008
    ),
    "At least 2 exposure variables are needed for variable selection"
  )
})

test_that("parallel Rho_G selection requires >=2 exposures even with multi-column CoG", {
  set.seed(1008)
  N <- 40
  G <- matrix(rnorm(N), nrow = N)            # one exposure only
  CoG <- matrix(rnorm(N * 2), nrow = N)      # extra covariates should not bypass rule
  Z1 <- matrix(rnorm(N * 4), nrow = N)
  Z2 <- matrix(rnorm(N * 4), nrow = N)
  Y <- rnorm(N)

  expect_error(
    estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      CoG = CoG,
      family = "normal",
      K = c(2, 2),
      Rho_G = 0.1,
      max_itr = 6,
      tol = 1e-1,
      seed = 1008
    ),
    "At least 2 exposure variables are needed for feature selection"
  )
})

test_that("parallel selectG stores overall union across layers", {
  set.seed(1008)
  G <- matrix(rnorm(200), nrow = 50)
  Z1 <- matrix(rnorm(300), nrow = 50)
  Z2 <- matrix(rnorm(300), nrow = 50)
  Y <- rnorm(50)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "parallel",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      family = "normal",
      K = c(2, 2),
      Rho_G = 0.1,
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  expect_type(fit$select$selectG, "logical")
  expect_true(is.list(fit$select$selectG_layer))
  expect_equal(
    fit$select$selectG,
    Reduce("|", fit$select$selectG_layer)
  )

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_parallel")
  expect_true(is.data.frame(s$feature_selection$G))
  expect_true(is.list(s$feature_selection$G_layer))
  expect_equal(
    s$feature_selection$G$Selected,
    Reduce("|", lapply(s$feature_selection$G_layer, `[[`, "Selected"))
  )
})

skip_on_cran()

# LUCID in serial - binary outcome (runtime-optimized coverage)

make_serial_tune_binary_data <- function(seed = 1008, n = 60, pG = 4, pZ = 6) {
  set.seed(seed)
  G <- matrix(rnorm(n * pG), nrow = n)
  Z1 <- matrix(rnorm(n * pZ), nrow = n)
  Z2 <- matrix(rnorm(n * pZ), nrow = n)
  Z3 <- matrix(rnorm(n * pZ), nrow = n)
  Z4 <- matrix(rnorm(n * pZ), nrow = n)
  Y <- rbinom(n = n, size = 1, prob = 0.30)
  CoY <- matrix(rnorm(n * 2), nrow = n)
  CoG <- matrix(rnorm(n * 2), nrow = n)
  list(G = G, Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Y = Y, CoY = CoY, CoG = CoG)
}

test_that("tune_lucid/lucid serial binary work for all-early topology", {
  d <- make_serial_tune_binary_data(seed = 1008)
  Z <- list(Z1 = d$Z1, Z2 = d$Z2, Z3 = d$Z3, Z4 = d$Z4)

  invisible(capture.output(
    tune <- tune_lucid(
      G = d$G, Z = Z, Y = d$Y,
      K = list(2:3, 2, 2, 2),
      lucid_model = "serial",
      family = "binary",
      init_omic.data.model = "VVV",
      CoG = d$CoG, CoY = d$CoY,
      seed = 1008,
      useY = TRUE
    )
  ))

  invisible(capture.output(
    fit <- lucid(
      G = d$G, Z = Z, Y = d$Y,
      K = list(2:3, 2, 2, 2),
      lucid_model = "serial",
      family = "binary",
      init_omic.data.model = "VVV",
      CoG = d$CoG, CoY = d$CoY,
      seed = 1008,
      useY = TRUE
    )
  ))

  expect_false(is.null(tune))
  expect_s3_class(fit, "lucid_serial")
})

test_that("tune_lucid/lucid serial binary work for mixed topology", {
  d <- make_serial_tune_binary_data(seed = 1010)
  Z <- list(Z1 = d$Z1, list(Z2 = d$Z2, Z3 = d$Z3), Z4 = d$Z4)

  invisible(capture.output(
    tune <- tune_lucid(
      G = d$G, Z = Z, Y = d$Y,
      K = list(2, list(2, 2:3), 2),
      lucid_model = "serial",
      family = "binary",
      init_omic.data.model = "VVV",
      CoG = d$CoG, CoY = d$CoY,
      seed = 1010,
      useY = TRUE
    )
  ))

  invisible(capture.output(
    fit <- lucid(
      G = d$G, Z = Z, Y = d$Y,
      K = list(2, list(2, 2:3), 2),
      lucid_model = "serial",
      family = "binary",
      init_omic.data.model = "VVV",
      CoG = d$CoG, CoY = d$CoY,
      seed = 1010,
      useY = TRUE
    )
  ))

  expect_false(is.null(tune))
  expect_s3_class(fit, "lucid_serial")
})

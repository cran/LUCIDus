# Serial stress smoke test: 6 consecutive parallel stages (2 layers each),
# with associated (non-independent) simulated data.

make_serial_six_parallel_data <- function(n = 48, pG = 5, pZ = 3, seed = 4242) {
  set.seed(seed)

  # Exposure block
  G <- matrix(rnorm(n * pG), nrow = n, ncol = pG)
  colnames(G) <- paste0("G", seq_len(pG))

  # Covariates derived from exposures to induce realistic association.
  CoG <- cbind(
    cov1 = G[, 1] + rnorm(n, sd = 0.25),
    cov2 = G[, 2] - 0.5 * G[, 3] + rnorm(n, sd = 0.25)
  )
  CoY <- CoG

  # Stage-specific latent signals with serial dependence.
  eta <- matrix(0, nrow = n, ncol = 6)
  x <- matrix(0, nrow = n, ncol = 6)

  eta[, 1] <- 0.9 * G[, 1] - 0.7 * G[, 2] + 0.4 * CoG[, 1] + rnorm(n, sd = 0.5)
  x[, 1] <- as.numeric(eta[, 1] > median(eta[, 1]))

  for (s in 2:6) {
    eta[, s] <- 0.7 * as.numeric(scale(eta[, s - 1])) +
      0.5 * G[, 1] - 0.4 * G[, 3] + 0.35 * x[, s - 1] + rnorm(n, sd = 0.6)
    x[, s] <- as.numeric(eta[, s] > median(eta[, s]))
  }

  # Build 6 serial stages, each stage is a 2-layer parallel block.
  Z <- vector("list", 6)
  names(Z) <- paste0("stage", seq_len(6))
  for (s in seq_len(6)) {
    shared <- rnorm(n, sd = 0.35)
    layer1 <- cbind(
      1.2 * x[, s] + 0.5 * G[, 1] + shared + rnorm(n, sd = 0.45),
      0.9 * x[, s] - 0.3 * G[, 2] + shared + rnorm(n, sd = 0.45),
      0.7 * x[, s] + 0.4 * G[, 4] + shared + rnorm(n, sd = 0.45)
    )
    layer2 <- cbind(
      -1.0 * x[, s] + 0.45 * G[, 2] + shared + rnorm(n, sd = 0.45),
      -0.8 * x[, s] - 0.35 * G[, 1] + shared + rnorm(n, sd = 0.45),
      -0.6 * x[, s] + 0.30 * G[, 5] + shared + rnorm(n, sd = 0.45)
    )
    colnames(layer1) <- paste0("s", s, "_L1_f", seq_len(pZ))
    colnames(layer2) <- paste0("s", s, "_L2_f", seq_len(pZ))
    Z[[s]] <- list(layer1 = layer1, layer2 = layer2)
  }

  # Outcome associated with final-stage latent signal + exposures + covariate.
  Y <- 0.8 * x[, 6] + 0.45 * G[, 1] - 0.25 * G[, 2] + 0.35 * CoY[, 2] +
    rnorm(n, sd = 0.7)

  list(G = G, Z = Z, Y = as.numeric(Y), CoG = CoG, CoY = CoY)
}

test_that("serial runs with 6 consecutive parallel submodels (K = 2 per layer)", {
  d <- make_serial_six_parallel_data(seed = 4242)
  K6 <- replicate(6, list(2, 2), simplify = FALSE)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = d$G,
      Z = d$Z,
      Y = d$Y,
      CoG = d$CoG,
      CoY = d$CoY,
      family = "normal",
      K = K6,
      Rho_G = 0,
      Rho_Z_Mu = 0,
      Rho_Z_Cov = 0,
      max_itr = 6,
      max_tot.itr = 30,
      tol = 1e-2,
      seed = 4242
    )
  )))

  expect_s3_class(fit, "lucid_serial")
  expect_equal(length(fit$submodel), 6)
  expect_true(all(vapply(fit$submodel, inherits, logical(1), what = "lucid_parallel")))
  expect_true(all(vapply(fit$submodel, function(sm) identical(as.numeric(sm$K), c(2, 2)), logical(1))))

  # Each parallel stage should provide 2-class inclusion probabilities per layer.
  expect_true(all(vapply(
    fit$submodel,
    function(sm) is.matrix(sm$inclusion.p[[1]]) && ncol(sm$inclusion.p[[1]]) == 2,
    logical(1)
  )))
  expect_true(all(vapply(
    fit$submodel,
    function(sm) is.matrix(sm$inclusion.p[[2]]) && ncol(sm$inclusion.p[[2]]) == 2,
    logical(1)
  )))

  s <- summary_lucid(fit)
  expect_s3_class(s, "sumlucid_serial")
  expect_equal(s$model_info$n_stages, 6)
  expect_equal(length(s$stage_summary), 6)
  expect_true(all(vapply(s$stage_summary, inherits, logical(1), what = "sumlucid_parallel")))
})

make_serial_six_parallel_monotone_missing_data <- function(n = 60, pG = 5, pZ = 3, seed = 5252) {
  d <- make_serial_six_parallel_data(n = n, pG = pG, pZ = pZ, seed = seed)

  # Monotone listwise missingness by stage:
  # once a row becomes all-missing in a stage, it stays all-missing in later stages.
  # Keep proportions moderate to avoid pathological initialization failures.
  listwise_n_by_stage <- c(0, 4, 8, 12, 16, 20)

  for (s in seq_len(6)) {
    if (listwise_n_by_stage[s] > 0) {
      miss_rows <- seq_len(listwise_n_by_stage[s])
      d$Z[[s]]$layer1[miss_rows, ] <- NA
      d$Z[[s]]$layer2[miss_rows, ] <- NA
    }

    # Add a small sporadic component that also increases with stage.
    spor_rows <- seq_len(min(n, s + 1))
    d$Z[[s]]$layer1[spor_rows, 1] <- NA
    d$Z[[s]]$layer2[spor_rows, pZ] <- NA
  }

  d$listwise_n_by_stage <- listwise_n_by_stage
  d
}

test_that("serial 6-stage parallel runs under monotone increasing stage-wise listwise missingness", {
  d <- make_serial_six_parallel_monotone_missing_data(seed = 5252)
  K6 <- replicate(6, list(2, 2), simplify = FALSE)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = d$G,
      Z = d$Z,
      Y = d$Y,
      CoG = d$CoG,
      CoY = d$CoY,
      family = "normal",
      K = K6,
      Rho_G = 0,
      Rho_Z_Mu = 0,
      Rho_Z_Cov = 0,
      max_itr = 6,
      max_tot.itr = 36,
      tol = 1e-2,
      seed = 5252
    )
  )))

  expect_s3_class(fit, "lucid_serial")
  expect_equal(length(fit$submodel), 6)
  expect_true(is.list(fit$missing_summary))
  expect_equal(fit$missing_summary$n_stages, 6)

  # Validate stage-wise monotone listwise pattern from recorded missing summaries.
  listwise_l1 <- vapply(fit$missing_summary$stage, function(ms) {
    as.integer(ms$layer_summary$listwise_rows[1])
  }, integer(1))
  listwise_l2 <- vapply(fit$missing_summary$stage, function(ms) {
    as.integer(ms$layer_summary$listwise_rows[2])
  }, integer(1))

  expect_equal(listwise_l1, d$listwise_n_by_stage)
  expect_equal(listwise_l2, d$listwise_n_by_stage)
  expect_true(all(diff(listwise_l1) >= 0))

  # Confirm model still produces finite probabilities in the final stage.
  last_stage <- fit$submodel[[6]]
  expect_true(all(is.finite(last_stage$inclusion.p[[1]])))
  expect_true(all(is.finite(last_stage$inclusion.p[[2]])))
})

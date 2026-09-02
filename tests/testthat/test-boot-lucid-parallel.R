skip_on_cran()

# Bootstrap smoke tests for parallel LUCID

test_that("boot_lucid parallel smoke test without covariates", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  expect_true(all(c("beta", "mu", "gamma", "bootstrap") %in% names(out)))
  expect_true(is.list(out$beta))
  expect_true(is.list(out$mu))
  expect_equal(length(out$beta), 2)
  expect_equal(length(out$mu), 2)
  expect_equal(nrow(out$beta[[1]]), (fit$K[1] - 1) * (ncol(G) + 1))
  expect_equal(nrow(out$beta[[2]]), (fit$K[2] - 1) * (ncol(G) + 1))
  expect_equal(nrow(out$mu[[1]]), fit$K[1] * ncol(Z1))
  expect_equal(nrow(out$mu[[2]]), fit$K[2] * ncol(Z2))
  expect_equal(ncol(out$beta[[1]]), 5)
  expect_equal(ncol(out$mu[[1]]), 5)

  gamma_len <- length(parallel_delta_coef(fit$res_Gamma$Gamma))
  expect_equal(nrow(out$gamma), gamma_len)
  expect_equal(ncol(out$gamma), 5)
  expect_true(is.finite(out$beta[[1]][1, "estimate"]))
  expect_true(is.finite(out$mu[[1]][1, "estimate"]))
  expect_s3_class(out$bootstrap, "boot")
})

test_that("boot_lucid parallel handles CoG and CoY indexing", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)
  CoG <- matrix(rnorm(40), nrow = 40)
  CoY <- matrix(rnorm(40), nrow = 40)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      CoG = CoG, CoY = CoY,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      CoG = CoG, CoY = CoY,
      lucid_model = "parallel",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  # per-layer beta block = intercept + exposures + CoG covariates
  n_beta_col <- ncol(G) + 1L + ncol(CoG)
  expect_equal(nrow(out$beta[[1]]), (fit$K[1] - 1) * n_beta_col)
  expect_equal(nrow(out$beta[[2]]), (fit$K[2] - 1) * n_beta_col)
  expect_equal(nrow(out$mu[[1]]), fit$K[1] * ncol(Z1))
  expect_equal(nrow(out$mu[[2]]), fit$K[2] * ncol(Z2))
})

test_that("boot_lucid parallel keeps exposure names with single G and CoG present", {
  set.seed(2026)
  G <- matrix(rnorm(40), nrow = 40)
  colnames(G) <- "hs_child_age_yrs_None"
  Z1 <- matrix(rnorm(120), nrow = 40)
  Z2 <- matrix(rnorm(120), nrow = 40)
  Y <- rnorm(40)
  CoG <- matrix(rnorm(40), nrow = 40)
  CoY <- matrix(rnorm(40), nrow = 40)
  colnames(CoG) <- "sex_male"
  colnames(CoY) <- "sex_male"

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      CoG = CoG, CoY = CoY,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      max_itr = 5,
      max_tot.itr = 80,
      tol = 2e-1,
      seed = 2026
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      CoG = CoG, CoY = CoY,
      lucid_model = "parallel",
      model = fit,
      R = 2,
      conf = 0.9
    )
  )))

  expect_true(any(grepl("hs_child_age_yrs_None", rownames(out$beta[[1]]), fixed = TRUE)))
  expect_true(any(grepl("hs_child_age_yrs_None", rownames(out$beta[[2]]), fixed = TRUE)))
  expect_true(is.finite(out$beta[[1]][1, "estimate"]))
  expect_true(is.finite(out$beta[[2]][1, "estimate"]))
})

test_that("boot_lucid parallel rejects models with unrefit feature selection", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  fit$select$selectG[1] <- FALSE

  expect_error(
    boot_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      model = fit,
      R = 2
    ),
    "Refit LUCID model with selected feature first"
  )
})

test_that("boot_lucid parallel integrates with summary_lucid and print", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    boot_out <- boot_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  s <- summary_lucid(fit, boot.se = boot_out)
  expect_s3_class(s, "sumlucid_parallel")
  expect_true(is.list(s$boot.se))
  expect_true(all(c("beta", "mu", "gamma", "bootstrap") %in% names(s$boot.se)))

  txt <- capture.output(print(s))
  expect_true(any(grepl("Detailed parameter estimates", txt)))
  expect_true(any(grepl("norm_lower", txt)))
  expect_true(any(grepl(rownames(boot_out$mu[[1]])[1], txt, fixed = TRUE)))
  beta_row_print <- sub("^Layer[0-9]+\\.", "", rownames(boot_out$beta[[1]])[1])
  expect_true(any(grepl(beta_row_print, txt, fixed = TRUE)))
})

test_that("parallel summary includes intercept in Y and E with and without bootstrap", {
  set.seed(1008)
  G <- matrix(rnorm(160), nrow = 40)
  colnames(G) <- paste0("g", 1:ncol(G))
  Z1 <- matrix(rnorm(320), nrow = 40)
  Z2 <- matrix(rnorm(320), nrow = 40)
  Y <- rnorm(40)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  s_plain <- summary_lucid(fit)
  txt_plain <- capture.output(print(s_plain))
  expect_true(any(grepl("^\\(Intercept\\)\\s", txt_plain)))
  expect_true(any(grepl("\\(Intercept\\)\\.cluster2", txt_plain)))

  suppressWarnings(invisible(capture.output(
    boot_out <- boot_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  s_boot <- summary_lucid(fit, boot.se = boot_out)
  txt_boot <- capture.output(print(s_boot))
  expect_true(any(grepl("^\\(Intercept\\)\\s", txt_boot)))
  expect_true(any(grepl("\\(Intercept\\)\\.cluster2", txt_boot)))
  expect_true(any(grepl("norm_lower", txt_boot)))
})

test_that("boot_lucid parallel runs with mixed missingness under mix imputation", {
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
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      init_impute = "mix",
      max_itr = 8,
      tol = 1e-1,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      model = fit,
      R = 2,
      conf = 0.9
    )
  )))

  expect_s3_class(out$bootstrap, "boot")
  expect_equal(length(out$beta), 2)
  expect_equal(length(out$mu), 2)
  expect_equal(ncol(out$gamma), 5)
})

test_that("boot_lucid parallel works for binary outcome", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000), nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Y <- rbinom(n = 100, size = 1, prob = 0.25)

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G,
      Z = list(Z1, Z2, Z3),
      Y = Y,
      lucid_model = "parallel",
      family = "binary",
      K = c(2, 2, 2),
      seed = i,
      useY = TRUE
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G,
      Z = list(Z1, Z2, Z3),
      Y = Y,
      lucid_model = "parallel",
      model = fit,
      R = 2,
      conf = 0.9
    )
  )))

  expect_true(all(c("beta", "mu", "gamma", "bootstrap") %in% names(out)))
  expect_equal(length(out$beta), 3)
  expect_equal(length(out$mu), 3)
  expect_equal(ncol(out$gamma), 5)
  expect_s3_class(out$bootstrap, "boot")
  expect_true(any(is.finite(out$gamma[, "estimate"])))
})

test_that("boot_lucid parallel auto-refits zero-penalty fallback when model has nonzero penalty", {
  set.seed(2027)
  G <- matrix(rnorm(180), nrow = 45)
  Z1 <- matrix(rnorm(225), nrow = 45)
  Z2 <- matrix(rnorm(225), nrow = 45)
  Y <- rnorm(45)

  suppressWarnings(invisible(capture.output(
    fit_pen <- estimate_lucid(
      G = G, Z = list(Z1, Z2), Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = c(2, 2),
      Rho_G = 0.01,
      Rho_Z_Mu = 0.01,
      Rho_Z_Cov = 0.01,
      max_itr = 8,
      max_tot.itr = 40,
      tol = 1e-1,
      seed = 2027
    )
  )))

  out <- NULL
  expect_warning(
    invisible(capture.output(
      out <- withCallingHandlers(
        boot_lucid(
          G = G, Z = list(Z1, Z2), Y = Y,
          lucid_model = "parallel",
          model = fit_pen,
          R = 2,
          conf = 0.9
        ),
        warning = function(w) {
          if (!grepl("zero-penalty", conditionMessage(w), fixed = TRUE)) {
            invokeRestart("muffleWarning")
          }
        }
      )
    )),
    "zero-penalty"
  )

  expect_true(all(c("beta", "mu", "gamma", "bootstrap") %in% names(out)))
  expect_equal(length(out$beta), 2)
  expect_equal(length(out$mu), 2)
  expect_true(is.matrix(out$gamma))
})

skip_on_cran()

# Bootstrap smoke tests for early LUCID

test_that("boot_lucid early smoke test without covariates", {
  G <- sim_data$G[1:60, ]
  Z <- sim_data$Z[1:60, ]
  Y <- sim_data$Y_normal[1:60, ]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  expect_true(all(c("beta", "mu", "gamma", "bootstrap") %in% names(out)))
  # beta block is the whole G->X matrix (intercept + exposures + CoG)
  expect_equal(nrow(out$beta), (fit$K - 1) * ncol(fit$res_Beta))
  expect_equal(nrow(out$mu), fit$K * ncol(Z))
  expect_equal(nrow(out$gamma), length(fit$res_Gamma$beta))
  expect_equal(ncol(out$beta), 5)
  expect_equal(ncol(out$mu), 5)
  expect_equal(ncol(out$gamma), 5)
  expect_s3_class(out$bootstrap, "boot")
})

test_that("boot_lucid early handles CoG and CoY indexing", {
  G <- sim_data$G[1:60, ]
  Z <- sim_data$Z[1:60, ]
  Y <- sim_data$Y_normal[1:60, ]
  CoG <- sim_data$Covariate[1:60, 1, drop = FALSE]
  CoY <- sim_data$Covariate[1:60, 2, drop = FALSE]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      CoG = CoG, CoY = CoY,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G, Z = Z, Y = Y,
      CoG = CoG, CoY = CoY,
      lucid_model = "early",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  # beta block now carries the intercept and CoG covariate rows too
  expect_equal(nrow(out$beta), (fit$K - 1) * ncol(fit$res_Beta))
  expect_equal(nrow(out$mu), fit$K * ncol(Z))
  expect_equal(nrow(out$gamma), length(fit$res_Gamma$beta))
})

test_that("boot_lucid early rejects models with unrefit feature selection", {
  G <- sim_data$G[1:60, ]
  Z <- sim_data$Z[1:60, ]
  Y <- sim_data$Y_normal[1:60, ]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  fit$select$selectG[1] <- FALSE

  expect_error(
    boot_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      model = fit,
      R = 2
    ),
    "Refit LUCID model with selected feature first"
  )
})

test_that("boot_lucid integrates with summary_lucid and print for early model", {
  G <- sim_data$G[1:60, ]
  Z <- sim_data$Z[1:60, ]
  Y <- sim_data$Y_normal[1:60, ]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    boot_out <- boot_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  s <- summary_lucid(fit, boot.se = boot_out)
  expect_s3_class(s, "sumlucid_early")
  expect_true(is.list(s$boot.se))
  expect_true(all(c("beta", "mu", "gamma", "bootstrap") %in% names(s$boot.se)))

  expect_no_error(capture.output(print(s)))
})

test_that("early summary uses consistent Y labels with and without bootstrap CIs", {
  G <- sim_data$G[1:80, ]
  Z <- sim_data$Z[1:80, ]
  Y <- sim_data$Y_normal[1:80, ]
  CoY <- sim_data$Covariate[1:80, , drop = FALSE]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      G = G, Z = Z, Y = Y, CoY = CoY,
      lucid_model = "early",
      family = "normal",
      K = 2,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    boot_out <- boot_lucid(
      G = G, Z = Z, Y = Y, CoY = CoY,
      lucid_model = "early",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  s_plain <- summary_lucid(fit)
  s_boot <- summary_lucid(fit, boot.se = boot_out)

  txt_plain <- capture.output(print(s_plain))
  txt_boot <- capture.output(print(s_boot))

  expect_true(any(grepl("^\\(Intercept\\)\\s", txt_plain)))
  expect_true(any(grepl("^\\(Intercept\\)\\s", txt_boot)))
  expect_true(any(grepl("^cluster2\\s", txt_plain)))
  expect_true(any(grepl("^cluster2\\s", txt_boot)))
  expect_true(any(grepl("^\\(Intercept\\)\\.cluster2\\s", txt_plain)))
  # the (3) E bootstrap CI table now shows the intercept row too
  expect_true(any(grepl("^\\(Intercept\\)\\.cluster2\\s", txt_boot)))
  expect_true(any(grepl(paste0("^", colnames(G)[1], "\\.cluster2\\s"), txt_boot)))
  expect_true(any(grepl("norm_lower", txt_boot)))
})

test_that("boot_lucid early auto-refits zero-penalty fallback when model has nonzero penalty", {
  G <- sim_data$G[1:70, ]
  Z <- sim_data$Z[1:70, ]
  Y <- sim_data$Y_normal[1:70, ]

  suppressWarnings(invisible(capture.output(
    fit_pen <- estimate_lucid(
      G = G, Z = Z, Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2,
      Rho_G = 0.01,
      Rho_Z_Mu = 0.01,
      Rho_Z_Cov = 0.01,
      max_itr = 8,
      max_tot.itr = 40,
      seed = 1016
    )
  )))

  out <- NULL
  expect_warning(
    invisible(capture.output(
      out <- withCallingHandlers(
        boot_lucid(
          G = G, Z = Z, Y = Y,
          lucid_model = "early",
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
  expect_equal(ncol(out$beta), 5)
  expect_equal(ncol(out$mu), 5)
  expect_equal(ncol(out$gamma), 5)
})

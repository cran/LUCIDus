# Bootstrap smoke tests for serial LUCID

test_that("boot_lucid serial smoke test for all-early topology", {
  G <- sim_data$G[1:50, ]
  Y <- sim_data$Y_normal[1:50, ]
  Z1 <- sim_data$Z[1:50, 1:4]
  Z2 <- sim_data$Z[1:50, 5:8]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = list(2, 2),
      family = "normal",
      max_itr = 4,
      max_tot.itr = 20,
      seed = 1008
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      lucid_model = "serial",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  expect_true(all(c("stage", "bootstrap") %in% names(out)))
  expect_s3_class(out$bootstrap, "boot")
  expect_equal(length(out$stage), 2)
  expect_true(all(c("beta", "mu", "gamma") %in% names(out$stage[[1]])))
  expect_true(all(c("beta", "mu", "gamma") %in% names(out$stage[[2]])))
  expect_true(is.null(out$stage[[1]]$gamma))
  expect_true(is.matrix(out$stage[[2]]$gamma))
  expect_equal(ncol(out$stage[[1]]$beta), 5)
  expect_equal(ncol(out$stage[[1]]$mu), 5)
})

test_that("boot_lucid serial integrates with summary_lucid and print for mixed topology", {
  G <- sim_data$G[1:50, ]
  Y <- sim_data$Y_normal[1:50, ]
  Z1 <- sim_data$Z[1:50, 1:3]
  Z2 <- sim_data$Z[1:50, 4:6]
  Z3 <- sim_data$Z[1:50, 7:10]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = list(list(Z1, Z2), Z3),
      Y = Y,
      K = list(list(2, 2), 2),
      family = "normal",
      max_itr = 4,
      max_tot.itr = 20,
      seed = 1010
    )
  )))

  suppressWarnings(invisible(capture.output(
    boot_out <- boot_lucid(
      G = G,
      Z = list(list(Z1, Z2), Z3),
      Y = Y,
      lucid_model = "serial",
      model = fit,
      R = 3,
      conf = 0.9
    )
  )))

  s <- summary_lucid(fit, boot.se = boot_out)
  expect_s3_class(s, "sumlucid_serial")
  expect_true(is.list(s$boot.se))
  expect_true(is.list(s$stage_summary[[1]]$boot.se))
  expect_true(is.list(s$stage_summary[[2]]$boot.se))

  txt <- capture.output(print(s))
  expect_equal(sum(grepl("^\\(1\\) Y \\(continuous outcome\\)", txt)), 1)
  expect_true(any(grepl("previous serial stage", txt, fixed = TRUE)))
  expect_true(any(grepl("norm_lower", txt, fixed = TRUE)))
})

test_that("boot_lucid serial rejects models with unrefit feature selection", {
  G <- sim_data$G[1:45, ]
  Y <- sim_data$Y_normal[1:45, ]
  Z1 <- sim_data$Z[1:45, 1:4]
  Z2 <- sim_data$Z[1:45, 5:8]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = list(2, 2),
      family = "normal",
      max_itr = 4,
      max_tot.itr = 20,
      seed = 1012
    )
  )))

  fit$submodel[[1]]$select$selectG[1] <- FALSE

  expect_error(
    boot_lucid(
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      lucid_model = "serial",
      model = fit,
      R = 2
    ),
    "Refit serial LUCID model with selected feature first"
  )
})

test_that("boot_lucid serial runs for binary outcome and keeps stage structure", {
  G <- sim_data$G[1:45, ]
  Y <- sim_data$Y_binary[1:45, ]
  Z1 <- sim_data$Z[1:45, 1:4]
  Z2 <- sim_data$Z[1:45, 5:8]

  suppressWarnings(invisible(capture.output(
    fit <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = list(2, 2),
      family = "binary",
      max_itr = 4,
      max_tot.itr = 20,
      seed = 1014
    )
  )))

  suppressWarnings(invisible(capture.output(
    out <- boot_lucid(
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      lucid_model = "serial",
      model = fit,
      R = 2,
      conf = 0.9
    )
  )))

  expect_equal(length(out$stage), 2)
  expect_true(is.null(out$stage[[1]]$gamma))
  expect_true(is.matrix(out$stage[[2]]$gamma))
  expect_equal(ncol(out$stage[[2]]$gamma), 5)
  expect_true(any(is.finite(out$stage[[2]]$gamma[, "estimate"])))
})

test_that("boot_lucid serial auto-refits zero-penalty fallback when model has nonzero penalty", {
  G <- sim_data$G[1:45, ]
  Y <- sim_data$Y_normal[1:45, ]
  Z1 <- sim_data$Z[1:45, 1:4]
  Z2 <- sim_data$Z[1:45, 5:8]

  suppressWarnings(invisible(capture.output(
    fit_pen <- estimate_lucid(
      lucid_model = "serial",
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      K = list(2, 2),
      family = "normal",
      Rho_G = 0.01,
      Rho_Z_Mu = 0.01,
      Rho_Z_Cov = 0.01,
      max_itr = 4,
      max_tot.itr = 20,
      seed = 1018
    )
  )))

  out <- NULL
  expect_warning(
    invisible(capture.output(
      out <- withCallingHandlers(
        boot_lucid(
          G = G,
          Z = list(Z1, Z2),
          Y = Y,
          lucid_model = "serial",
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

  expect_true(all(c("stage", "bootstrap") %in% names(out)))
  expect_equal(length(out$stage), 2)
  expect_true(is.null(out$stage[[1]]$gamma))
  expect_true(is.matrix(out$stage[[2]]$gamma))
})

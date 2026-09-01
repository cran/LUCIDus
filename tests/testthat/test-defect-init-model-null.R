# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# init_omic.data.model = NULL is the documented way to let mclust choose the
# covariance geometry, and is what the R Journal vignette passes. For the
# parallel model it used to fail: R/em.R built the per-layer model names with
# rep(init_omic.data.model, length(K)), and rep(NULL, n) is NULL, so the names
# were lost between initialization and the M-step. mclust::mstep() then received
# modelName = NULL and died in switch(EXPR = modelName, ...) with "EXPR must be
# a length 1 vector", which reached the user as the far less informative
# "LUCID model fails to converge given current tuning parameters".

test_that("automatic covariance selection works for the parallel model", {
  for (n_layer in c(1L, 2L, 3L)) {
    d <- sim_lucid("parallel", N = 300, K = 2, M = 4, n_layer = n_layer, seed = 401)
    suppressWarnings(suppressMessages(invisible(capture.output(
      fit <- estimate_lucid(lucid_model = "parallel", G = d$G, Z = d$Z, Y = d$Y,
                            family = "normal", K = rep(2, n_layer),
                            init_omic.data.model = NULL)
    ))))
    expect_s3_class(fit, "lucid_parallel")
    expect_sane_fit(fit, info = paste("parallel NULL", n_layer, "layer"))
    # the selected geometry is reported back, one name per layer, rather than
    # left as the NULL that was asked for
    expect_length(fit$init_omic.data.model, n_layer)
    expect_true(all(nzchar(fit$init_omic.data.model)))
  }
})

test_that("automatic covariance selection works through the lucid() wrapper", {
  d <- sim_lucid("parallel", N = 300, K = 2, M = 4, n_layer = 2, seed = 402)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "parallel",
                 family = "normal", K = list(2, 2), init_omic.data.model = NULL)
  ))))
  expect_s3_class(fit, "lucid_parallel")
  expect_sane_fit(fit, info = "parallel NULL via lucid()")
})

test_that("automatic covariance selection works for serial, including a parallel stage", {
  topo <- list(2, list(2, 2))
  d <- sim_lucid_serial_topology(topo, N = 300, M = 4, seed = 403)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial",
                 family = "normal", K = topo, init_omic.data.model = NULL)
  ))))
  expect_s3_class(fit, "lucid_serial")
  expect_sane_fit(fit, info = "serial NULL with parallel stage")
})

test_that("a single geometry name applies to every layer", {
  d <- sim_lucid("parallel", N = 300, K = 2, M = 4, n_layer = 3, seed = 404)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(lucid_model = "parallel", G = d$G, Z = d$Z, Y = d$Y,
                          family = "normal", K = c(2, 2, 2),
                          init_omic.data.model = "EII")
  ))))
  expect_equal(fit$init_omic.data.model, rep("EII", 3))
})

test_that("a per-layer geometry vector is recycled to one name per layer", {
  # The same line used rep() where rep_len() was meant, so a length-n vector for
  # n layers came back with length n^2. It happened to index correctly, so
  # nothing failed -- but the wrong-length vector was stored on the fit.
  d <- sim_lucid("parallel", N = 300, K = 2, M = 4, n_layer = 2, seed = 405)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(lucid_model = "parallel", G = d$G, Z = d$Z, Y = d$Y,
                          family = "normal", K = c(2, 2),
                          init_omic.data.model = c("EII", "VVV"))
  ))))
  expect_equal(fit$init_omic.data.model, c("EII", "VVV"))
})

test_that("random initialization falls back to a concrete geometry", {
  # Random initialization has no mclust fit to read a selection from, so NULL
  # must resolve to the same EEV default the early path uses rather than
  # reaching the M-step unset.
  d <- sim_lucid("parallel", N = 300, K = 2, M = 4, n_layer = 2, seed = 406)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- estimate_lucid(lucid_model = "parallel", G = d$G, Z = d$Z, Y = d$Y,
                          family = "normal", K = c(2, 2),
                          init_omic.data.model = NULL, init_par = "random")
  ))))
  expect_equal(fit$init_omic.data.model, rep("EEV", 2))
  expect_sane_fit(fit, info = "parallel NULL random init")
})

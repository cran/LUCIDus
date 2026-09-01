skip_on_cran()

# Combination sweep: every architecture x topology x family x missingness cell
# fits, returns the structure it was asked for, and contains no NA, NaN, Inf or
# absurd parameter value.
#
# This is the breadth tier. It does not check accuracy -- test-accuracy-*.R does
# that -- but it is the guarantee that no supported combination is simply
# broken, which is exactly the class of defect that reached the vignettes twice.
#
# Cells are exhaustive on the axes that interact in the source (architecture x
# layer count, serial topology, family, missingness) and fixed at defaults
# elsewhere.

N_GRID <- 300L
M_GRID <- 4L

fit_cell <- function(model, K, d, family, ...) {
  suppressWarnings(suppressMessages(invisible(capture.output(
    fit <- lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = model,
                 family = family, K = K, init_omic.data.model = NULL, ...)
  ))))
  fit
}

# ---- early -----------------------------------------------------------------

test_that("early fits every family x missingness cell and returns sane values", {
  for (family in c("normal", "binary")) {
    for (miss in c("none", "listwise", "sporadic", "mixed")) {
      lbl <- paste("early", family, miss)
      d <- sim_lucid("early", N = N_GRID, K = 2, M = M_GRID,
                     family = family, missing = miss, seed = 101)
      fit <- fit_cell("early", 2, d, family)

      expect_s3_class(fit, "early_lucid")
      expect_sane_fit(fit, info = lbl)
      expect_ascends(fit, info=lbl)
      expect_equal(nrow(fit$inclusion.p), N_GRID)
      expect_equal(ncol(fit$inclusion.p), 2L)
      expect_equal(dim(fit$res_Mu), c(2L, M_GRID))
      if (miss != "none") expect_listwise_preserved(fit$Z, d$miss_index, info = lbl)
    }
  }
})

# ---- parallel: 1, 2, 3 and 5 layers ----------------------------------------

test_that("parallel fits 1, 2, 3 and 5 layers for both families", {
  for (family in c("normal", "binary")) {
    for (n_layer in c(1L, 2L, 3L, 5L)) {
      lbl <- paste("parallel", family, n_layer, "layer")
      d <- sim_lucid("parallel", N = N_GRID, K = 2, M = M_GRID,
                     n_layer = n_layer, family = family, seed = 102)
      fit <- fit_cell("parallel", as.list(rep(2, n_layer)), d, family)

      expect_s3_class(fit, "lucid_parallel")
      expect_sane_fit(fit, info = lbl)
      expect_ascends(fit, info=lbl)
      # one posterior block per layer, each N x K
      expect_length(fit$inclusion.p, n_layer)
      for (m in fit$inclusion.p) expect_equal(dim(as.matrix(m)), c(N_GRID, 2L))
      expect_length(fit$res_Mu, n_layer)
      expect_equal(fit$N, N_GRID)
    }
  }
})

test_that("parallel tolerates every missingness regime at 2 and 3 layers", {
  for (n_layer in c(2L, 3L)) {
    for (miss in c("listwise", "sporadic", "mixed")) {
      lbl <- paste("parallel", n_layer, "layer", miss)
      d <- sim_lucid("parallel", N = N_GRID, K = 2, M = M_GRID,
                     n_layer = n_layer, missing = miss, seed = 103)
      fit <- fit_cell("parallel", as.list(rep(2, n_layer)), d, "normal")

      expect_sane_fit(fit, info = lbl)

      expect_ascends(fit, info=lbl)
      expect_listwise_preserved(fit$Z, d$miss_index, info = lbl)
    }
  }
})

# ---- serial: stage counts --------------------------------------------------

test_that("serial fits 1, 2, 3 and 6 stages for both families", {
  for (family in c("normal", "binary")) {
    for (n_stage in c(1L, 2L, 3L, 6L)) {
      lbl <- paste("serial", family, n_stage, "stage")
      d <- sim_lucid("serial", N = N_GRID, K = 2, M = M_GRID,
                     n_layer = n_stage, family = family, seed = 104)
      Z <- if (n_stage == 1L) list(d$Z[[1]]) else d$Z
      fit <- fit_cell("serial", as.list(rep(2, n_stage)),
                      list(G = d$G, Z = Z, Y = d$Y), family)

      expect_s3_class(fit, "lucid_serial")
      expect_sane_fit(fit, info = lbl)
      expect_ascends(fit, info=lbl)
      expect_length(fit$submodel, n_stage)
      # one transition per consecutive pair, so none at all for a single stage
      expect_length(fit$res_Delta, max(0L, n_stage - 1L))
    }
  }
})

# ---- serial: stage topology, the "early and parallel within" requirement ----

test_that("serial fits every mix of early and parallel stages", {
  topologies <- list(
    "all early"            = list(2, 2),
    "all parallel"         = list(list(2, 2), list(2, 2)),
    "early then parallel"  = list(2, list(2, 2)),
    "parallel then early"  = list(list(2, 2), 2),
    "early parallel early" = list(2, list(2, 2), 2),
    "three stage mixed"    = list(list(2, 2), 2, list(2, 2))
  )
  for (nm in names(topologies)) {
    topo <- topologies[[nm]]
    d <- sim_lucid_serial_topology(topo, N = N_GRID, M = M_GRID, seed = 105)
    fit <- fit_cell("serial", topo, d, "normal")

    expect_s3_class(fit, "lucid_serial")
    expect_sane_fit(fit, info = nm)
    expect_ascends(fit, info=nm)
    expect_length(fit$submodel, length(topo))
    # a stage given a list of cluster counts must be fitted as a parallel
    # sub-model, and one given a single count as an early sub-model
    expected_class <- vapply(topo, function(k) {
      if (is.list(k)) "lucid_parallel" else "early_lucid"
    }, character(1))
    expect_equal(vapply(fit$submodel, function(s) class(s)[1], character(1)),
                 expected_class, info = nm)
  }
})

test_that("serial mixed topology survives missing data and both families", {
  topo <- list(2, list(2, 2), 2)
  for (family in c("normal", "binary")) {
    for (miss in c("listwise", "sporadic", "mixed")) {
      lbl <- paste("serial mixed", family, miss)
      d <- sim_lucid_serial_topology(topo, N = N_GRID, M = M_GRID,
                                     family = family, missing = miss, seed = 106)
      fit <- fit_cell("serial", topo, d, family)
      expect_sane_fit(fit, info = lbl)
      expect_ascends(fit, info=lbl)
      expect_length(fit$submodel, 3L)
    }
  }
})

# ---- K > 2, unequal K, covariates, supervision ------------------------------

test_that("unequal cluster counts across layers and stages are honoured", {
  d <- sim_lucid("parallel", N = N_GRID, K = c(2, 3), M = M_GRID,
                 n_layer = 2, seed = 107)
  fit <- fit_cell("parallel", list(2, 3), d, "normal")
  expect_sane_fit(fit, info = "parallel K=2,3")
  expect_ascends(fit, info="parallel K=2,3")
  expect_equal(ncol(as.matrix(fit$inclusion.p[[1]])), 2L)
  expect_equal(ncol(as.matrix(fit$inclusion.p[[2]])), 3L)

  d2 <- sim_lucid_serial_topology(list(3, list(2, 2)), N = N_GRID, M = M_GRID, seed = 108)
  fit2 <- fit_cell("serial", list(3, list(2, 2)), d2, "normal")
  expect_sane_fit(fit2, info = "serial K=3 then 2,2")
  expect_ascends(fit2, info="serial K=3 then 2,2")
  expect_equal(ncol(as.matrix(fit2$submodel[[1]]$inclusion.p)), 3L)
})

test_that("covariate adjustment on either arm fits for every architecture", {
  # CoG adjusts G -> X, CoY adjusts X -> Y; they enter through different code
  # paths and are worth crossing with architecture.
  d_e <- sim_lucid("early", N = N_GRID, K = 2, M = M_GRID, n_CoG = 2, n_CoY = 2, seed = 109)
  fit_e <- suppressWarnings(suppressMessages(invisible(capture.output(
    f <- lucid(G = d_e$G, Z = d_e$Z, Y = d_e$Y, CoG = d_e$CoG, CoY = d_e$CoY,
               lucid_model = "early", family = "normal", K = 2,
               init_omic.data.model = NULL)))))
  expect_sane_fit(f, info = "early CoG+CoY")
  expect_ascends(f, info="early CoG+CoY")

  d_p <- sim_lucid("parallel", N = N_GRID, K = 2, M = M_GRID, n_layer = 2,
                   n_CoG = 2, n_CoY = 2, seed = 110)
  suppressWarnings(suppressMessages(invisible(capture.output(
    fp <- lucid(G = d_p$G, Z = d_p$Z, Y = d_p$Y, CoG = d_p$CoG, CoY = d_p$CoY,
                lucid_model = "parallel", family = "normal", K = list(2, 2),
                init_omic.data.model = NULL)))))
  expect_sane_fit(fp, info = "parallel CoG+CoY")
  expect_ascends(fp, info="parallel CoG+CoY")
})

test_that("unsupervised fits are sane for every architecture", {
  for (spec in list(list("early", 2, 1L), list("parallel", list(2, 2), 2L),
                    list("serial", list(2, 2), 2L))) {
    model <- spec[[1]]; K <- spec[[2]]; n_layer <- spec[[3]]
    d <- sim_lucid(model, N = N_GRID, K = 2, M = M_GRID, n_layer = n_layer, seed = 111)
    suppressWarnings(suppressMessages(invisible(capture.output(
      fit <- estimate_lucid(lucid_model = model, G = d$G, Z = d$Z, Y = d$Y,
                            family = "normal",
                            K = if (model == "early") 2 else rep(2, n_layer),
                            useY = FALSE, init_omic.data.model = NULL)))))
    expect_sane_fit(fit, info = paste(model, "unsupervised"))
    expect_ascends(fit, info=paste(model, "unsupervised"))
    expect_false(fit$useY)
  }
})

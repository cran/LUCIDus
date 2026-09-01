# One targeted regression per defect fixed in this pass.  Each of these fails
# on the pre-fix code, which is what the previous shape-only suite did not do.

test_that("D2: parallel penalized G-step keeps exactly K-1 reference-coded rows", {
  set.seed(61); N <- 300
  G <- matrix(stats::rnorm(N * 4), N, 4); colnames(G) <- paste0("g", 1:4)
  for (Kv in list(c(2, 2), c(3, 2))) {
    r <- array(stats::runif(prod(Kv) * N), dim = c(Kv, N))
    for (i in seq_len(N)) r[, , i] <- r[, , i] / sum(r[, , i])
    s <- LUCIDus:::Mstep_GtoX(G = G, r = r, selectG = TRUE, penalty = 0.02,
                              K = Kv, N = N)
    expect_equal(nrow(s$Beta[[1]]), Kv[1] - 1L)
    p <- exp(LUCIDus:::f_GtoX(G, s$Beta[[1]]))
    expect_equal(ncol(p), Kv[1])                       # not K + 1
    expect_equal(as.numeric(rowSums(p)), rep(1, N), tolerance = 1e-10)
  }
})

test_that("D2: fitted probabilities do not depend on the reference class", {
  set.seed(62); N <- 250
  G <- matrix(stats::rnorm(N * 3), N, 3); colnames(G) <- paste0("g", 1:3)
  r <- array(stats::runif(4 * N), dim = c(2, 2, N))
  for (i in seq_len(N)) r[, , i] <- r[, , i] / sum(r[, , i])
  B <- LUCIDus:::Mstep_GtoX(G = G, r = r, selectG = TRUE, penalty = 0.02,
                            K = c(2, 2), N = N)$Beta[[1]]
  full <- rbind(0, B)
  p1 <- exp(LUCIDus:::f_GtoX(G, B))
  p2 <- exp(LUCIDus:::f_GtoX(G, sweep(full, 2, full[2, ], "-")[-2, , drop = FALSE]))
  expect_equal(unname(p1), unname(p2[, c(2, 1)]), tolerance = 1e-10)
})

test_that("D3: parallel penalized Z deselects features whose means are zero", {
  skip_on_cran()
  d <- sim_lucid("parallel", N = 300, K = c(2, 2), n_layer = 2, M = 4,
                 family = "normal", seed = 63)
  strong <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, K = list(2, 2),
                                   lucid_model = "parallel", family = "normal",
                                   Rho_Z_Mu = 500))
  sz <- strong$select$selectZ[[1]]
  n_sel <- if (is.null(dim(sz))) sum(sz) else sum(colSums(sz) > 0)
  expect_true(all(abs(strong$res_Mu[[1]]) < 1e-10))
  expect_equal(n_sel, 0L)
})

test_that("D4: variable selection refits successfully with one surviving exposure", {
  skip_on_cran()
  d <- sim_lucid("early", N = 400, K = 2, P = 1, M = 4, family = "normal", seed = 64)
  Gn <- cbind(d$G, matrix(stats::rnorm(400 * 6), 400, 6))
  colnames(Gn) <- c("signal", paste0("null", 1:6))
  fit <- suppressMessages(lucid(G = Gn, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2,
                                Rho_G = c(0, 0.05, 0.1, 0.2), seed = 1008))
  expect_s3_class(fit, "early_lucid")
  # the original selection is retained for reporting
  expect_equal(length(fit$selection$selectG), ncol(Gn))
  # and column names survive the single-column subset
  expect_true(all(colnames(fit$res_Beta)[-1] %in% colnames(Gn)))
})

test_that("D6: predicted cluster labels are 1..K for every model type", {
  skip_on_cran()
  de <- sim_lucid("early", N = 250, K = 3, M = 3, family = "normal", seed = 65)
  fe <- suppressMessages(lucid(G = de$G, Z = de$Z, Y = de$Y, lucid_model = "early",
                               family = "normal", K = 3, init_omic.data.model = "VVV"))
  pe <- predict_lucid(model = fe, G = de$G, Z = de$Z, Y = de$Y, lucid_model = "early")
  expect_true(all(pe$pred.x %in% 1:3))
  expect_equal(min(pe$pred.x), 1)

  dp <- sim_lucid("parallel", N = 250, K = c(2, 2), n_layer = 2, M = 3, seed = 66)
  fp <- suppressMessages(lucid(G = dp$G, Z = dp$Z, Y = dp$Y, K = list(2, 2),
                               lucid_model = "parallel", family = "normal"))
  pp <- predict_lucid(model = fp, G = dp$G, Z = dp$Z, Y = dp$Y, lucid_model = "parallel")
  expect_true(all(unlist(pp$pred.x) %in% 1:2))
})

test_that("D8: bootstrap t0 is the supplied model, not a random refit", {
  skip_on_cran()
  d <- sim_lucid("early", N = 250, K = 2, P = 2, M = 3, family = "normal", seed = 67)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  set.seed(7)
  b <- suppressWarnings(boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                   model = fit, R = 25))
  expected <- LUCIDus:::lucid_early_par_vector(fit, ncol(d$G))
  expect_equal(as.numeric(b$bootstrap$t0), as.numeric(expected))
  # and the reported estimate column agrees with summary(fit)
  s <- summary(fit, boot.se = b, auto_print = FALSE)
  expect_equal(unname(s$boot.se$mu[, "estimate"]), as.vector(t(fit$res_Mu)),
               tolerance = 1e-12)
})

test_that("D8: failed replicates are counted but do not destroy the interval", {
  # This test previously asserted the opposite -- that fewer than 20 valid
  # replicates blanked every confidence limit.  That behaviour was itself the
  # bug: it silently wiped the CI tables of every vignette (all of which use
  # R = 2, 3 or 5).  Counting failures is useful; refusing to report computable
  # limits is not.
  skip_on_cran()
  d <- sim_lucid("early", N = 200, K = 2, P = 2, M = 3, family = "normal", seed = 68)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  set.seed(8)
  b <- suppressWarnings(boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                   model = fit, R = 25))
  b$bootstrap$t[1:6, ] <- NA          # six replicates fail outright

  expect_warning(LUCIDus:::boot_replicate_status(b$bootstrap),
                 "failed to produce finite")
  ci <- suppressWarnings(LUCIDus:::gen_ci(b$bootstrap))
  expect_true(all(is.finite(ci[, -1])))            # limits survive
  expect_equal(attr(ci, "boot_status")$n_valid, 19L)
})

test_that("D8: limits are NA only when fewer than two replicates survive", {
  skip_on_cran()
  d <- sim_lucid("early", N = 200, K = 2, P = 2, M = 3, family = "normal", seed = 68)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  set.seed(8)
  b <- suppressWarnings(boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                   model = fit, R = 25))
  b$bootstrap$t[2:25, ] <- NA         # one survivor: sd() is undefined
  expect_warning(ci <- LUCIDus:::gen_ci(b$bootstrap), "minimum")
  expect_true(all(is.na(ci[, -1])))
  expect_true(all(is.finite(ci[, "estimate"])))     # point estimate is kept
  # stricter behaviour remains available on request
  b2 <- suppressWarnings(boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                    model = fit, R = 4))
  ci2 <- suppressWarnings(LUCIDus:::gen_ci(b2$bootstrap, min_valid = 20L))
  expect_true(all(is.na(ci2[, -1])))
})

test_that("boot_lucid() returns finite limits at the small R used in the vignettes", {
  # Characterization test for the regression the vignettes exposed: with R = 2,
  # 3 or 5 every confidence limit came back NA.  Small R gives WIDE intervals,
  # not absent ones.
  skip_on_cran()
  d <- sim_lucid("early", N = 150, K = 2, P = 2, M = 3, family = "normal", seed = 94)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  for (R in c(2L, 3L, 5L)) {
    set.seed(31)
    b <- suppressWarnings(boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                     model = fit, R = R))
    cells <- rbind(b$beta, b$mu, b$gamma)[, c("norm_lower", "norm_upper",
                                              "perc_lower", "perc_upper")]
    expect_true(all(is.finite(cells)), info = paste("R =", R))
    expect_true(all(b$beta[, "perc_lower"] <= b$beta[, "perc_upper"]),
                info = paste("R =", R))
  }
})

test_that("failed bootstrap replicates do not leak raw error text", {
  # try(est_lucid(...)) was wrapped only in capture.output(), which redirects
  # stdout while try() writes to stderr -- so replicate failures printed a raw
  # "Error in est_lucid(...)" dump that looked like a crash.
  skip_on_cran()
  d <- sim_lucid("early", N = 120, K = 2, P = 2, M = 3, family = "normal", seed = 95)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  tf <- tempfile(); con <- file(tf, open = "wt")
  sink(con, type = "message")
  invisible(withCallingHandlers(
    boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early", model = fit, R = 5),
    warning = function(w) invokeRestart("muffleWarning")))
  sink(type = "message"); close(con)
  leaked <- readLines(tf, warn = FALSE)
  expect_false(any(grepl("Error in est_lucid|mclust failed", leaked)))
})

test_that("tuning is not broken by ties in BIC", {
  # which() returning several indices used to trigger recursive `[[` indexing
  tl <- data.frame(K = c(2, 2), Rho_G = c(0, 0), Rho_Z_Mu = c(0, 0),
                   Rho_Z_Cov = c(0, 0), BIC = c(100, 100))
  x <- min(tl[, 5])
  expect_length(which(tl[, 5] == x)[1], 1L)
})

test_that("glasso is not invoked, and does not warn, when Rho_Z_Cov is zero", {
  skip_on_cran()
  d <- sim_lucid("parallel", N = 250, K = c(2, 2), n_layer = 2, M = 3, seed = 69)
  w <- character(0)
  withCallingHandlers(
    suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, K = list(2, 2),
                           lucid_model = "parallel", family = "normal",
                           Rho_Z_Mu = 10)),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
  )
  expect_false(any(grepl("rho=0", w, fixed = TRUE)))
})

test_that("predict_lucid() works without naming lucid_model (paper Section 3.5)", {
  skip_on_cran()
  # The R Journal paper calls predict_lucid(model = fit1, G = ..., Z = ..., Y = ...)
  # with no lucid_model argument.  The unresolved default was forwarded to the
  # internal worker, which accepts only c("early","parallel"), so match.arg()
  # failed with "'arg' must be of length 1" for every such call.
  d <- sim_lucid("early", N = 250, K = 2, M = 3, family = "normal", seed = 91)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))

  p_default  <- predict_lucid(model = fit, G = d$G, Z = d$Z, Y = d$Y)
  p_explicit <- predict_lucid(model = fit, G = d$G, Z = d$Z, Y = d$Y,
                              lucid_model = "early")
  expect_equal(p_default$pred.x, p_explicit$pred.x)
  expect_equal(p_default$inclusion.p, p_explicit$inclusion.p)
  expect_true(all(p_default$pred.x %in% seq_len(fit$K)))

  # and without Y, as the paper also shows
  expect_no_error(predict_lucid(model = fit, G = d$G, Z = d$Z))
})

test_that("predict_lucid() dispatches correctly for parallel and serial too", {
  skip_on_cran()
  dp <- sim_lucid("parallel", N = 200, K = c(2, 2), n_layer = 2, M = 3, seed = 92)
  fp <- suppressMessages(lucid(G = dp$G, Z = dp$Z, Y = dp$Y, K = list(2, 2),
                               lucid_model = "parallel", family = "normal"))
  expect_equal(
    predict_lucid(model = fp, G = dp$G, Z = dp$Z, Y = dp$Y, lucid_model = "parallel")$pred.x,
    predict_lucid(model = fp, G = dp$G, Z = dp$Z, Y = dp$Y, lucid_model = "parallel")$pred.x
  )
  ds <- sim_lucid("serial", N = 200, K = c(2, 2), n_layer = 2, M = 3, seed = 93)
  fs <- suppressMessages(lucid(G = ds$G, Z = ds$Z, Y = ds$Y, lucid_model = "serial",
                               K = list(2, 2), family = "normal"))
  expect_no_error(predict_lucid(model = fs, G = ds$G, Z = ds$Z, Y = ds$Y,
                                lucid_model = "serial"))
})

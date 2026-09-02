# Heavy: fits LUCID models + bootstraps; runs locally and in CI, not on CRAN.
skip_on_cran()

# Regression tests for the 2026 correctness pass:
#  - boot_lucid() aligns each replicate's cluster labels to the reference fit,
#    so a coefficient whose canonical ordering statistic is weak no longer gets
#    a bimodal / sign-split bootstrap interval (early / parallel / serial).
#  - summary() bootstrap CI tables carry an odds-ratio view and headers that
#    match the printed scale.
#  - "Finished ... log-likelihood" completion messages, and loglik_trace, report
#    the observed-data (unpenalized) value.

# ------------------------------------------------------------------ helpers ----
mk_binary <- function(n = 240, sep_col = 1, layers = 1, seed = 11) {
  set.seed(seed)
  xt <- rbinom(n, 1, 0.5)
  G <- cbind(g1 = (xt - 0.5) * 1.3 + rnorm(n, sd = 0.7), g2 = rnorm(n))
  mkZ <- function() {
    Z <- matrix(rnorm(n * 4), n, 4)
    Z[, sep_col] <- Z[, sep_col] + 3 * xt          # strong omics separation
    colnames(Z) <- paste0("z", 1:4)
    Z
  }
  Y <- as.integer((0.2 + 1.5 * xt + rnorm(n)) > 0)
  Z <- if (layers == 1) mkZ() else replicate(layers, mkZ(), simplify = FALSE)
  list(G = G, Z = Z, Y = matrix(Y, ncol = 1), xt = xt)
}

t0_inside_ci <- function(b, conf = 0.95) {
  tt <- b$bootstrap$t
  t0 <- b$bootstrap$t0
  a <- (1 - conf) / 2
  ok <- rep(TRUE, length(t0))
  for (j in seq_along(t0)) {
    v <- tt[, j]; v <- v[is.finite(v)]
    if (length(v) < 3) next
    lo <- stats::quantile(v, a); hi <- stats::quantile(v, 1 - a)
    ok[j] <- t0[j] >= lo - 1e-6 && t0[j] <= hi + 1e-6
  }
  ok
}

sign_stability <- function(b, col) {
  tt <- b$bootstrap$t
  nm <- names(b$bootstrap$t0)
  v <- tt[, which(nm == col)]; v <- v[is.finite(v)]
  max(mean(v > 0), mean(v < 0))
}

# ------------------------------------------------------------- parallel --------
test_that("parallel: every replicate coefficient's t0 lies inside its own CI", {
  d <- mk_binary(n = 220, layers = 2, seed = 21)
  fit <- suppressWarnings(suppressMessages(estimate_lucid(
    G = d$G, Z = d$Z, Y = d$Y, lucid_model = "parallel", family = "binary",
    K = c(2, 2), seed = 1, max_itr = 60, max_tot.itr = 120, tol = 1e-3)))
  set.seed(99)
  b <- suppressWarnings(suppressMessages(
    boot_lucid(G = d$G, Z = d$Z, Y = d$Y, model = fit, R = 60, conf = 0.95)))

  expect_true(all(t0_inside_ci(b)),
              info = "aligned replicates keep the point estimate inside its bootstrap interval")
  # the strong G->X coefficient must be one-signed across replicates
  expect_gt(sign_stability(b, "Layer1.g1.cluster2"), 0.9)
})

# --------------------------------------------------------------- serial --------
test_that("serial: transition and outcome coefficients bracket t0 after alignment", {
  d <- mk_binary(n = 200, layers = 2, seed = 31)
  fit <- suppressWarnings(suppressMessages(estimate_lucid(
    G = d$G, Z = d$Z, Y = d$Y, lucid_model = "serial", family = "binary",
    K = list(2, 2), seed = 1, max_itr = 40, max_tot.itr = 80, tol = 1e-3)))
  set.seed(99)
  b <- suppressWarnings(suppressMessages(
    boot_lucid(G = d$G, Z = d$Z, Y = d$Y, model = fit, R = 40, conf = 0.95)))

  expect_true(all(t0_inside_ci(b)))
  # the final-stage outcome log-OR is one-signed
  nm <- names(b$bootstrap$t0)
  y_col <- nm[grepl("Y\\.LC2$", nm)][1]
  expect_false(is.na(y_col))
  expect_gt(sign_stability(b, y_col), 0.9)
})

# ---------------------------------------------------------------- early --------
test_that("early: alignment leaves the outcome-anchored fit and t0 untouched", {
  d <- mk_binary(n = 220, layers = 1, seed = 41)
  fit <- suppressWarnings(suppressMessages(estimate_lucid(
    G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early", family = "binary",
    K = 2, seed = 1, max_itr = 60, max_tot.itr = 120, tol = 1e-3)))
  set.seed(7)
  b <- suppressWarnings(suppressMessages(
    boot_lucid(G = d$G, Z = d$Z, Y = d$Y, model = fit, R = 40, conf = 0.95)))

  # Part 1c: the reference side (t0) is never permuted
  expect_equal(as.numeric(b$bootstrap$t0),
               as.numeric(LUCIDus:::lucid_early_par_vector(fit, ncol(d$G))))
  expect_true(all(t0_inside_ci(b)))
})

test_that("boot_lucid is deterministic given the seed, with or without lucid_model", {
  d <- mk_binary(n = 180, layers = 1, seed = 51)
  fit <- suppressWarnings(suppressMessages(estimate_lucid(
    G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early", family = "binary",
    K = 2, seed = 1, max_itr = 40, max_tot.itr = 80, tol = 1e-3)))
  set.seed(3)
  b1 <- suppressWarnings(suppressMessages(
    boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early", model = fit, R = 15)))
  set.seed(3)
  b2 <- suppressWarnings(suppressMessages(
    boot_lucid(G = d$G, Z = d$Z, Y = d$Y, model = fit, R = 15)))
  expect_equal(b1$bootstrap$t, b2$bootstrap$t)
})

# --------------------------------------------------- CI-table OR + wording -----
test_that("early binary CI table gains an OR view and a scale-honest header", {
  d <- mk_binary(n = 220, layers = 1, seed = 61)
  fit <- suppressWarnings(suppressMessages(estimate_lucid(
    G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early", family = "binary",
    K = 2, seed = 1, max_itr = 60, max_tot.itr = 120, tol = 1e-3)))
  set.seed(5)
  b <- suppressWarnings(suppressMessages(
    boot_lucid(G = d$G, Z = d$Z, Y = d$Y, model = fit, R = 40, conf = 0.95)))

  txt_plain <- utils::capture.output(print(summary(fit, auto_print = FALSE)))
  txt_boot  <- utils::capture.output(print(summary(fit, boot.se = b, auto_print = FALSE)))

  # (3) E header: no-CI branch promises an odds ratio (an OR column is shown);
  # CI branch says "coefficient" scale and mentions the interval.
  e_plain <- grep("^\\(3\\) E:", txt_plain, value = TRUE)
  e_boot  <- grep("^\\(3\\) E:", txt_boot,  value = TRUE)
  expect_match(e_plain, "odds ratio")
  expect_match(e_boot, "coefficient")
  expect_match(e_boot, "interval")

  # binary (1) Y CI table carries OR / OR_lower / OR_upper
  expect_true(any(grepl("\\bOR_lower\\b", txt_boot)))
  expect_true(any(grepl("\\bOR_upper\\b", txt_boot)))
})

test_that(".with_or_columns exponentiates the coefficient-scale columns", {
  df <- data.frame(estimate = c(0, log(2)), norm_lower = c(-1, 0.1),
                   norm_upper = c(1, 1.2), sig = c("", "*"),
                   check.names = FALSE)
  out <- LUCIDus:::.with_or_columns(df)
  expect_equal(out$OR, exp(df$estimate))
  expect_equal(out$OR_lower, exp(df$norm_lower))
  expect_equal(out$OR_upper, exp(df$norm_upper))
  expect_identical(names(out)[length(names(out))], "sig")   # sig stays last
})

# ------------------------------------------------- log-likelihood labeling -----
test_that("penalized fit reports the observed-data log-likelihood", {
  d <- mk_binary(n = 200, layers = 1, seed = 71)
  msg <- utils::capture.output(
    fit <- suppressWarnings(estimate_lucid(
      G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early", family = "binary",
      K = 2, Rho_G = 0.05, seed = 1, verbose = TRUE,
      max_itr = 40, max_tot.itr = 80, tol = 1e-2)))
  finished <- grep("^Finished LUCID early model:", msg, value = TRUE)
  expect_match(finished, "observed-data log-likelihood")
  expect_false(any(grepl("Finished.*penalized log-likelihood", msg)))

  # trace is the unpenalized observed loglik -> its last value is $likelihood
  tr <- fit$em_control$loglik_trace
  expect_gt(length(tr), 0)
  expect_equal(tr[length(tr)], fit$likelihood, tolerance = 1e-6)
})

# ---------------------------------------- (3) E intercept + covariate CI rows ---
# Part 2a: the bootstrap (3) E table must show a CI for the intercept and every
# CoG covariate, so it has the same coefficient rows as summary(fit) without a
# bootstrap. Checked for all three model types and both outcome families; the
# (1) Y intercept + CoY rows (already working) are a regression guard.

mk_cov_data <- function(n = 150, layers = 1, family = "normal", seed = 81) {
  set.seed(seed)
  G <- matrix(rnorm(n * 4), n, 4); colnames(G) <- paste0("g", 1:4)
  Z <- matrix(rnorm(n * 6), n, 6); colnames(Z) <- paste0("z", 1:6)
  CoG <- matrix(rnorm(n), n, 1); colnames(CoG) <- "age"
  CoY <- matrix(rnorm(n), n, 1); colnames(CoY) <- "bmi"
  lp <- 0.7 * G[, 1] - 0.6 * G[, 2] + rnorm(n)
  Y <- if (family == "binary") matrix(as.integer(lp > median(lp)), ncol = 1) else matrix(lp, ncol = 1)
  Zl <- if (layers == 1) Z else lapply(seq_len(layers), function(j) Z[, 1:3] + j)
  list(G = G, Z = Zl, Y = Y, CoG = CoG, CoY = CoY)
}

# distinct E-block coefficient row tokens ("<feature>.cluster<k>", with the
# optional "Layer<i>." prefix stripped) across every "(N) E" block. The E
# features in mk_cov_data() are (Intercept), g1..g4, age, and the serial
# transition Stage1.cluster<k>; (2) Z rows (z1..z6) are excluded.
e_tokens <- function(txt) {
  m <- regmatches(txt, regexpr("^([A-Za-z(][^ ]*?)\\.cluster[0-9]+", txt))
  m <- sub("^Layer[0-9]+\\.", "", m[nzchar(m)])
  keep <- grepl("^(\\(Intercept\\)|g[0-9]+|age|Stage[0-9]+\\.cluster[0-9]+)\\.cluster[0-9]+$", m)
  sort(unique(m[keep]))
}

for (mt in c("early", "parallel", "serial")) {
  for (fam in c("normal", "binary")) {
    test_that(sprintf("(3) E bootstrap CI shows intercept + CoG rows [%s / %s]", mt, fam), {
      d <- mk_cov_data(layers = if (mt == "early") 1 else 2, family = fam,
                       seed = 80 + which(c("early","parallel","serial") == mt) * 3 +
                              which(c("normal","binary") == fam))
      K <- switch(mt, early = 2, parallel = c(2, 2), serial = list(2, 2))
      fit <- suppressWarnings(suppressMessages(estimate_lucid(
        G = d$G, Z = d$Z, Y = d$Y, CoG = d$CoG, CoY = d$CoY,
        lucid_model = mt, family = fam, K = K, seed = 1,
        max_itr = 25, max_tot.itr = 60, tol = 1e-2)))
      set.seed(4)
      b <- suppressWarnings(suppressMessages(boot_lucid(
        G = d$G, Z = d$Z, Y = d$Y, CoG = d$CoG, CoY = d$CoY,
        lucid_model = mt, model = fit, R = 30, conf = 0.9)))

      txt_plain <- utils::capture.output(print(summary(fit, auto_print = FALSE)))
      txt_boot  <- utils::capture.output(print(summary(fit, boot.se = b, auto_print = FALSE)))

      tok_boot  <- e_tokens(txt_boot)
      tok_plain <- e_tokens(txt_plain)
      expect_true(any(grepl("^\\(Intercept\\)\\.cluster", tok_boot)),
                  info = "(3) E CI table has an intercept row")
      expect_true(any(grepl("^age\\.cluster", tok_boot)),
                  info = "(3) E CI table has the CoG covariate row")
      # exactly the same coefficient rows with and without the bootstrap CI
      expect_setequal(tok_boot, tok_plain)

      # (1) Y intercept + CoY still present (regression guard)
      y_start <- grep("^\\([0-9]+\\) Y", txt_boot)
      expect_true(length(y_start) > 0)
      y_blk <- txt_boot[seq(y_start[length(y_start)], length(txt_boot))]
      expect_true(any(grepl("Intercept", y_blk)))
      expect_true(any(grepl("bmi", y_blk)))
    })
  }
}

# --------------------------- serial + binary (1) Y bootstrap CI is not all-NA ---
# f.binary.early() used to join the CI to x$beta by row name; a serial model's
# final stage prefixes the outcome coefficient names with "Y.", so every match()
# was NA and the whole (1) Y CI table (limits, OR bounds, sig) came out NA.
test_that("serial binary (1) Y bootstrap CI table has finite limits", {
  set.seed(11); n <- 170
  xt <- rbinom(n, 1, 0.5)
  G  <- scale(cbind(g1 = (xt - .5) * 1.4 + rnorm(n, sd = .6), g2 = rnorm(n)))
  Z1 <- matrix(rnorm(n * 4), n, 4); Z1[, 1:2] <- Z1[, 1:2] + 3 * xt; colnames(Z1) <- paste0("a", 1:4)
  Z2 <- matrix(rnorm(n * 4), n, 4); Z2[, 1:2] <- Z2[, 1:2] + 3 * xt; colnames(Z2) <- paste0("b", 1:4)
  CoY <- matrix(rnorm(n), n, 1); colnames(CoY) <- "bmi"
  Yb <- as.integer((0.3 + 1.5 * xt + rnorm(n)) > 0)

  fit <- suppressWarnings(suppressMessages(estimate_lucid(
    G = G, Z = list(Z1, Z2), Y = matrix(Yb, ncol = 1), CoY = CoY,
    lucid_model = "serial", family = "binary", K = list(2, 2),
    seed = 1, max_itr = 25, max_tot.itr = 60, tol = 1e-2)))
  set.seed(7)
  b <- suppressWarnings(suppressMessages(boot_lucid(
    G = G, Z = list(Z1, Z2), Y = matrix(Yb, ncol = 1), CoY = CoY,
    lucid_model = "serial", model = fit, R = 30, conf = 0.9)))

  txt <- utils::capture.output(print(summary(fit, boot.se = b, auto_print = FALSE)))
  y <- grep("^\\(1\\) Y \\(binary outcome\\)", txt)
  expect_true(length(y) > 0)
  blk <- txt[seq(y[length(y)], min(y[length(y)] + 12, length(txt)))]
  # the (Intercept) / LC2 / bmi rows must carry finite numeric CI bounds, not NA
  rows <- grep("^(\\(Intercept\\)|LC2|bmi)\\s", blk, value = TRUE)
  expect_gte(length(rows), 3)
  nums <- unlist(regmatches(rows, gregexpr("-?[0-9]+\\.[0-9]+(e[+-][0-9]+)?", rows)))
  expect_true(length(nums) >= length(rows) * 3)          # gamma + at least 2 CI cols each
  expect_false(any(grepl("\\bNA\\b", rows)))
})

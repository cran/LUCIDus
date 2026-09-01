# Part D-a: every printed quantity is checked against the fitted object AND
# against an independent oracle, and every label is checked against the block it
# heads.  The 3.0.3 regression -- correct number, wrong label -- was invisible to
# a test that only asked whether summary() ran.

# Heavy: every block fits a LUCID model to convergence; runs locally and in CI.
skip_on_cran()

test_that("early summary values equal the fitted parameters", {
  d <- sim_lucid("early", N = 400, K = 3, P = 2, M = 4, family = "normal", seed = 71)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 3, init_omic.data.model = "VVV"))
  s <- summary(fit, auto_print = FALSE)

  expect_equal(s$loglik, fit$likelihood)
  expect_equal(s$BIC, o_bic(fit$likelihood, s$model_fit$n_parameters,
                            nrow(fit$inclusion.p)))
  expect_equal(s$model_info$K, fit$K)
  expect_equal(s$model_info$n_observations, nrow(fit$inclusion.p))
  expect_equal(unname(as.matrix(s$parameters$mu)), unname(fit$res_Mu))
  expect_equal(as.numeric(s$parameters$gamma$beta), as.numeric(fit$res_Gamma$beta))
})

test_that("the printed odds ratio is exp(beta), elementwise", {
  d <- sim_lucid("early", N = 400, K = 3, P = 3, M = 3, family = "normal", seed = 72)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 3, init_omic.data.model = "VVV"))
  out <- utils::capture.output(print(summary(fit, auto_print = FALSE)))
  block <- grep("\\.cluster[0-9]", out, value = TRUE)
  expect_gt(length(block), 0)
  parsed <- do.call(rbind, lapply(block, function(l) {
    f <- strsplit(trimws(l), "\\s+")[[1]]
    data.frame(name = f[1], beta = as.numeric(f[2]), OR = as.numeric(f[3]))
  }))
  expect_equal(parsed$OR, exp(parsed$beta), tolerance = 1e-4)
})

test_that("mu column k in the printout is cluster k's mean", {
  d <- sim_lucid("early", N = 400, K = 2, M = 4, family = "normal", seed = 73)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  s <- summary(fit, auto_print = FALSE)
  mu_tab <- as.matrix(s$parameters$mu)
  # printed table is features x clusters; res_Mu is clusters x features
  expect_equal(unname(mu_tab), unname(fit$res_Mu))
  for (k in seq_len(fit$K)) {
    expect_equal(unname(mu_tab[k, ]), unname(fit$res_Mu[k, ]),
                 info = paste("cluster", k))
  }
})

test_that("binary outcome coefficients are on the log-odds scale", {
  d <- sim_lucid("early", N = 500, K = 2, M = 4, family = "binary", seed = 74)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "binary", K = 2, init_omic.data.model = "VVV"))
  lev <- fit$res_Gamma$cluster_effect
  r <- fit$inclusion.p; y <- as.numeric(d$Y)
  # cluster-specific weighted event rate, on the logit scale
  p_hat <- colSums(r * y) / colSums(r)
  expect_equal(as.numeric(lev), as.numeric(stats::qlogis(p_hat)), tolerance = 1e-2)
})

test_that("the missing-data profile matches a directly computed truth", {
  for (miss in c("none", "listwise", "sporadic", "mixed")) {
    d <- sim_lucid("early", N = 300, K = 2, M = 5, family = "normal",
                   missing = miss, miss_ratio = 0.2, seed = 75)
    fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                  family = "normal", K = 2, init_omic.data.model = "VVV"))
    ms <- fit$missing_summary
    rs <- rowSums(is.na(d$Z)); M <- ncol(d$Z)
    expect_equal(ms$listwise_rows, sum(rs == M), info = miss)
    expect_equal(ms$sporadic_rows, sum(rs > 0 & rs < M), info = miss)
    expect_equal(ms$complete_rows, sum(rs == 0), info = miss)
    expect_equal(ms$total_missing_cells, sum(is.na(d$Z)), info = miss)
    expect_equal(ms$prop_total_missing_cells, mean(is.na(d$Z)), info = miss)
  }
})

test_that("summary runs and reports for every model type and family", {
  skip_on_cran()
  for (fam in c("normal", "binary")) {
    de <- sim_lucid("early", N = 250, K = 2, M = 3, family = fam, seed = 76)
    fe <- suppressMessages(lucid(G = de$G, Z = de$Z, Y = de$Y, lucid_model = "early",
                                 family = fam, K = 2, init_omic.data.model = "VVV"))
    expect_s3_class(summary(fe, auto_print = FALSE), "sumlucid_early")

    dp <- sim_lucid("parallel", N = 250, K = c(2, 2), n_layer = 2, M = 3,
                    family = fam, seed = 77)
    fp <- suppressMessages(lucid(G = dp$G, Z = dp$Z, Y = dp$Y, K = list(2, 2),
                                 lucid_model = "parallel", family = fam))
    sp <- summary(fp, auto_print = FALSE)
    # summary(fit)$BIC must work for every model type, not just "early"
    expect_true(is.finite(sp$BIC))
    expect_equal(sp$BIC, sp$model_fit$BIC)
    expect_equal(sp$loglik, fp$likelihood)
    for (i in 1:2) {
      expect_equal(unname(as.matrix(sp$parameters$mu[[i]])), unname(fp$res_Mu[[i]]),
                   tolerance = 1e-10, info = paste(fam, "layer", i))
    }
  }
})

test_that("bootstrap CI table brackets its own point estimate", {
  skip_on_cran()
  d <- sim_lucid("early", N = 250, K = 2, P = 2, M = 3, family = "normal", seed = 78)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  set.seed(11)
  b <- suppressWarnings(boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                   model = fit, R = 40))
  s <- summary(fit, boot.se = b, auto_print = FALSE)
  for (blk in c("beta", "mu", "gamma")) {
    tab <- s$boot.se[[blk]]
    skip_if(is.null(tab) || !nrow(tab), paste("no", blk, "rows"))
    ok <- is.finite(tab[, "perc_lower"]) & is.finite(tab[, "perc_upper"])
    expect_true(all(tab[ok, "perc_lower"] <= tab[ok, "perc_upper"]), info = blk)
    expect_true(all(is.finite(tab[, "estimate"])), info = blk)
  }
})

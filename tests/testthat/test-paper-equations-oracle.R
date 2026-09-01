# Heavy: fits multiple LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# Tier 1: analytic checks of every estimating equation against independent
# oracles in helper-oracle.R.  No EM run to convergence, so these are fast.

test_that("Eq 2: reported log-likelihood equals the observed-data likelihood", {
  d <- sim_lucid("early", N = 300, K = 2, P = 2, M = 4, family = "normal", seed = 21)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  orc <- o_estep_early(d$G, d$Z, d$Y, fit$res_Beta, fit$res_Mu,
                       simplify2array(fit$res_Sigma),
                       fit$res_Gamma$cluster_effect, fit$res_Gamma$sigma,
                       "normal", TRUE)
  expect_equal(fit$likelihood, orc$loglik, tolerance = 1e-6)
})

test_that("Eqs 3/4: responsibilities equal the oracle posterior", {
  d <- sim_lucid("early", N = 300, K = 2, M = 4, family = "normal", seed = 22)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  orc <- o_estep_early(d$G, d$Z, d$Y, fit$res_Beta, fit$res_Mu,
                       simplify2array(fit$res_Sigma),
                       fit$res_Gamma$cluster_effect, fit$res_Gamma$sigma,
                       "normal", TRUE)
  expect_equal(unname(fit$inclusion.p), unname(orc$r), tolerance = 1e-6)
  expect_posterior_valid(fit$inclusion.p)
})

test_that("Eqs 6/7: mu and Sigma are responsibility-weighted moments", {
  d <- sim_lucid("early", N = 400, K = 2, M = 4, family = "normal", seed = 23)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  r <- fit$inclusion.p
  mu_oracle <- (t(r) %*% d$Z) / colSums(r)
  expect_equal(unname(fit$res_Mu), unname(mu_oracle), tolerance = 1e-3)

  # VVV is the unconstrained model, so Eq 7 holds exactly
  for (k in 1:2) {
    S <- crossprod(sweep(d$Z, 2, fit$res_Mu[k, ]) * sqrt(r[, k])) / sum(r[, k])
    expect_equal(unname(fit$res_Sigma[[k]]), unname(S), tolerance = 1e-2)
  }
})

test_that("Eqs 8/10: gamma and sigma are responsibility-weighted outcome moments", {
  d <- sim_lucid("early", N = 400, K = 2, M = 4, family = "normal", seed = 24)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  r <- fit$inclusion.p; y <- as.numeric(d$Y)
  lev <- colSums(r * y) / colSums(r)
  sdv <- sqrt(colSums(r * sapply(lev, function(m) (y - m)^2)) / colSums(r))
  expect_equal(as.numeric(fit$res_Gamma$cluster_effect), as.numeric(lev), tolerance = 1e-3)
  expect_equal(as.numeric(fit$res_Gamma$sigma), as.numeric(sdv), tolerance = 1e-3)
})

test_that("Eq 12: unsupervised responsibilities ignore Y", {
  d <- sim_lucid("early", N = 300, K = 2, M = 4, family = "normal", seed = 25)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, useY = FALSE,
                                init_omic.data.model = "VVV"))
  orc <- o_estep_early(d$G, d$Z, d$Y, fit$res_Beta, fit$res_Mu,
                       simplify2array(fit$res_Sigma), NULL, NULL,
                       "normal", useY = FALSE)
  expect_equal(unname(fit$inclusion.p), unname(orc$r), tolerance = 1e-6)
})

test_that("Eq 13: unpenalized BIC matches the parameter count formula", {
  d <- sim_lucid("early", N = 300, K = 2, P = 2, M = 4, family = "normal", seed = 26)
  for (fam in c("normal", "binary")) {
    Y <- if (fam == "normal") d$Y else matrix(as.numeric(d$Y > mean(d$Y)), ncol = 1)
    for (K in c(2, 3)) {
      fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = Y, lucid_model = "early",
                                    family = fam, K = K, init_omic.data.model = "VVV"))
      s <- summary(fit, auto_print = FALSE)
      expect_equal(s$model_fit$n_parameters,
                   o_npars_early(ncol(d$G), ncol(d$Z), K, fam),
                   info = paste(fam, "K =", K))
      expect_equal(s$BIC, o_bic(fit$likelihood, s$model_fit$n_parameters, nrow(d$G)))
    }
  }
})

test_that("Eq 18: penalized BIC subtracts one per deselected variable", {
  d <- sim_lucid("early", N = 400, K = 2, P = 2, M = 4, family = "normal", seed = 27)
  Gn <- cbind(d$G, matrix(stats::rnorm(400 * 4), 400, 4))
  colnames(Gn) <- c(colnames(d$G), paste0("null", 1:4))
  fit <- suppressMessages(lucid(G = Gn, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, Rho_G = 0.05, seed = 5))
  s <- summary(fit, auto_print = FALSE)
  DG <- sum(!fit$select$selectG); DZ <- sum(!fit$select$selectZ)
  expect_equal(s$model_fit$n_parameters,
               o_npars_early(length(fit$select$selectG), ncol(fit$Z), fit$K,
                             "normal", DG = DG, DZ = DZ))
})

test_that("Eq 21: predicted cluster is the argmax of the posterior, on 1..K", {
  d <- sim_lucid("early", N = 300, K = 3, M = 4, family = "normal", seed = 28)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 3, init_omic.data.model = "VVV"))
  p <- predict_lucid(model = fit, G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early")
  expect_equal(p$pred.x, as.numeric(max.col(p$inclusion.p)))
  expect_true(all(p$pred.x %in% seq_len(fit$K)))
  expect_posterior_valid(p$inclusion.p)
})

test_that("Eq 17: the I-step closed form matches the oracle and is a maximizer", {
  set.seed(31)
  M <- 4; K <- 2
  mu <- sim_default_mu(K, M, sep = 1.2)
  sigma <- array(0, dim = c(M, M, K))
  for (k in 1:K) {
    A <- matrix(stats::rnorm(M * M), M)
    sigma[, , k] <- crossprod(A) / M + diag(0.5, M)
  }
  obs <- as.numeric(sim_rmvnorm(1, mu[1, ], sigma[, , 1]))
  obs[c(2, 4)] <- NA
  r <- c(0.3, 0.7)

  pkg <- fill_data(obs = obs, mu = t(mu), sigma = sigma, p = r,
                   index = !is.na(obs), lucid_model = "early")
  orc <- o_istep_row(obs, mu, sigma, r)
  expect_equal(as.numeric(pkg), as.numeric(orc), tolerance = 1e-8)

  # Eq 17 is a fixed-point (majorization) step: it maximizes the weighted
  # Gaussian log-density with the mixing weights HELD FIXED at their current
  # values w_k = r_k * phi(Z^(t) | mu_k, Sigma_k).  It is not the stationary
  # point of the self-consistent objective in which the weights also move with
  # z, so the check must hold the weights fixed to be a fair one.
  # p_ij in Eq 17 is phi(Z_i^(t) | mu_j, Sigma_j): the FULL density at the
  # current iterate, whose missing entries are the responsibility-weighted mean
  # (the value fill_data() starts from) -- not the marginal on observed
  # coordinates.  Using the marginal here would define a different objective.
  cur0 <- obs
  cur0[is.na(obs)] <- as.numeric(r %*% mu[, is.na(obs), drop = FALSE])
  dens0 <- vapply(1:K, function(k) exp(o_dmvnorm_obs(cur0, mu[k, ], sigma[, , k])), numeric(1))
  w_fixed <- r * dens0
  w_fixed <- w_fixed / sum(w_fixed)
  qfun <- function(z) {
    full <- obs; full[is.na(obs)] <- z
    sum(w_fixed * vapply(1:K, function(k) o_dmvnorm_obs(full, mu[k, ], sigma[, , k]),
                         numeric(1)))
  }
  opt <- stats::optim(c(0, 0), function(z) -qfun(z), method = "BFGS",
                      control = list(reltol = 1e-14))
  expect_gte(qfun(pkg[is.na(obs)]) + 1e-6, -opt$value)
  # and the closed form should BE the optimizer, not merely as good
  expect_equal(as.numeric(pkg[is.na(obs)]), as.numeric(opt$par), tolerance = 1e-4)
})

test_that("Eqs 19/20: bootstrap CIs match hand-computed intervals", {
  skip_on_cran()
  d <- sim_lucid("early", N = 250, K = 2, P = 2, M = 3, family = "normal", seed = 32)
  fit <- suppressMessages(lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                family = "normal", K = 2, init_omic.data.model = "VVV"))
  set.seed(99)
  b <- suppressWarnings(boot_lucid(G = d$G, Z = d$Z, Y = d$Y, lucid_model = "early",
                                   model = fit, R = 40))
  ci <- LUCIDus:::gen_ci(b$bootstrap, conf = 0.95, min_valid = 5L)
  for (i in seq_len(min(4, length(b$bootstrap$t0)))) {
    hand <- o_boot_ci(b$bootstrap$t0[i], b$bootstrap$t[, i], conf = 0.95)
    expect_equal(unname(ci[i, c("norm_lower", "norm_upper")]),
                 unname(hand$normal), tolerance = 1e-6)
    expect_equal(unname(ci[i, c("perc_lower", "perc_upper")]),
                 unname(hand$percentile), tolerance = 1e-6)
  }
})

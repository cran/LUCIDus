# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# LUCID - 1 omics, binary outcome

test_that("check estimations of LUCID with binary outcome (K = 2)", {


  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  i <- 1008

  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))
  invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = 2,
                                             lucid_model = "early",
                                             family = "binary",
                                             CoY = cov,

                                             seed = i,
                                             useY = TRUE)))

  summary(fit1)
  pars <- fit1
  beta_causal <- mean(abs(pars$res_Beta[2, 2:5]))
  beta_non <- mean(abs(pars$res_Beta[2, 6:10]))
  mu_causal <- mean(abs(pars$res_Mu[1, 1:5] - pars$res_Mu[2, 1:5]))
  mu_non <- mean(abs(pars$res_Mu[1, 6:10] - pars$res_Mu[2, 6:10]))
  gamma <- as.numeric(pars$res_Gamma$beta)

  # check parameters via signal-vs-noise separation (more stable than exact values)
  expect_true(is.finite(beta_causal))
  expect_true(is.finite(beta_non))
  expect_gt(beta_causal, beta_non)
  expect_gt(mu_causal, mu_non)
  expect_equal(length(gamma), 4)
  expect_true(all(is.finite(gamma)))
  expect_true(sum(abs(gamma[-1]) > 1e-6) >= 1)

  # check summary_lucid
  sum_fit1 <- summary(fit1)
  expect_equal(class(fit1), "early_lucid")
  expect_equal(class(sum_fit1), "sumlucid_early")

  # check predict_lucid early pred.x contract
  pred <- predict_lucid(
    model = fit1,
    lucid_model = "early",
    G = G,
    Z = Z,
    Y = Y,
    CoY = cov,
    response = FALSE
  )
  expect_type(pred$pred.x, "double")
  expect_equal(length(pred$pred.x), nrow(G))
  # D6: cluster labels are 1..K, matching Eq 21 and every other output block
  expect_true(all(pred$pred.x %in% seq_len(fit1$K)))
  expect_true(all(is.finite(pred$pred.x)))
  expect_true(is.matrix(pred$inclusion.p))
  expect_equal(dim(pred$inclusion.p), c(nrow(G), fit1$K))
  expect_equal(rowSums(pred$inclusion.p), rep(1, nrow(G)), tolerance = 1e-6)
})

test_that("check variable selection on G", {
  # run LUCID model
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  # cat(paste("test2 - seed =", i, "\n"))
  i <- 1008
  suppressWarnings(invisible(capture.output(fit1 <- lucid(G = G,
                                         Z = Z,
                                         Y = Y_binary,
                                         CoY = cov,
                                         family = "binary",
                                         lucid_model = "early",
                                         K = 2,
                                         seed = i,
                                         useY = TRUE,

                                         Rho_G = 0.1))))

  # check parameters
  expect_equal(class(fit1$select$selectG), "logical")
  expect_equal(as.vector(fit1$select$selectG),
               rep(TRUE, 4))
})


test_that("check variable selection on Z", {
  # run LUCID model
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  # cat(paste("test3 - seed =", i, "\n"))
  i <- 1008
  suppressWarnings(invisible(capture.output(fit1 <- lucid(G = G,
                                         Z = Z,
                                         Y = Y_binary,
                                         CoY = cov,
                                         family = "binary",
                                         lucid_model = "early",
                                         K = 2,
                                         seed = i,
                                         useY = TRUE,

                                         Rho_Z_Mu =  50,
                                         Rho_Z_Cov = 0.5))))

  # check parameters
  expect_equal(class(fit1$select$selectZ), "logical")
})


test_that("check whether arguments of lucid work", {
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  # cat(paste("test4 - seed =", i, "\n"))
  i <- 1008
  suppressWarnings(invisible(capture.output(fit1 <- lucid(G = G,
                                         Z = Z,
                                         Y = Y_binary,
                                         CoY = cov,
                                         family = "binary",
                                         lucid_model = "early",
                                         K = 2,
                                         seed = i,
                                         useY = TRUE,
                                         init_omic.data.model = NULL))))
  suppressWarnings(invisible(capture.output(fit2 <- lucid(G = G,
                                         Z = Z,
                                         Y = Y_binary,
                                         CoY = cov,
                                         family = "binary",
                                         lucid_model = "early",
                                         K = 2,
                                         seed = i,
                                         useY = TRUE,
                                         init_omic.data.model = "EEV",
                                         init_par = "random"))))
  expect_equal(class(fit1$init_omic.data.model), "character")
  expect_equal(fit2$init_par, "random")
})


test_that("check whether lucid throws an erorr with continuous outcome or outcome coded rather than 0 and 1", {
  G <- sim_data$G[1:500, ]
  Z <- sim_data$Z[1:500, ]
  Y_binary <- sim_data$Y_binary[1:500, ]
  Y_normal <- sim_data$Y_normal[1:500, ]
  cov <- sim_data$Covariate[1:500, ]
  # i <- sample(1:2000, 1)
  i <- 1008
  # cat(paste("test4 - seed =", i, "\n"))
  expect_error(est_lucid(G = G,
                         Z = Z,
                         Y = Y_normal,
                         CoY = cov,
                         family = "binary",
                         lucid_model = "early",
                         K = 2,
                         seed = i,
                         useY = TRUE,
                         init_omic.data.model = NULL))
  expect_error(est_lucid(G = G,
                         Z = Z,
                         Y = Y_binary + 1,
                         CoY = cov,
                         family = "binary",
                         lucid_model = "early",
                         K = 2,
                         seed = i,
                         useY = TRUE,
                         init_omic.data.model = NULL))
})

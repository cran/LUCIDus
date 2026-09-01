# Heavy: fits LUCID models; runs locally and in CI, not on CRAN.
skip_on_cran()

# LUCID - LUCID in Serial, binary outcome



test_that("check estimations of LUCID in Serial with binary outcome (K = 2,2,2,2,2)", {
  i <- 1008
  set.seed(i)
  G <- matrix(rnorm(500), nrow = 100)
  Z1 <- matrix(rnorm(1000),nrow = 100)
  Z2 <- matrix(rnorm(1000), nrow = 100)
  Z3 <- matrix(rnorm(1000), nrow = 100)
  Z4 <- matrix(rnorm(1000), nrow = 100)
  Z5 <- matrix(rnorm(1000), nrow = 100)

  ##missing data
  a = sample(1:1000, 30, replace=FALSE)
  Z1[a] = NA
  Z2[62:65, 6:8] = NA
  a = sample(1:1000, 30, replace=FALSE)
  Z4[a] = NA

  Z <- list(list(Z1 = Z1, Z2 = Z2, Z3 = Z3), Z4 = Z4, Z5 = Z5)
  Y <- rbinom(n=100, size =1, prob =0.25)
  CoY <- matrix(rnorm(200), nrow = 100)
  CoG <- matrix(rnorm(200), nrow = 100)

  # i <- sample(1:2000, 1)
  # cat(paste("test1 - seed =", i, "\n"))

  ##needs work!!!!!#####
  ##needs work!!!!!#####
  ##needs work!!!!!#####
  numeric_leaf <- function(x) {
    if (is.numeric(x)) return(as.numeric(x))
    if (is.list(x)) return(unlist(lapply(x, numeric_leaf), use.names = FALSE))
    numeric(0)
  }

  suppressWarnings(invisible(capture.output(fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(list(2,2,2),2,2),
                                             lucid_model = "serial",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             CoG = CoG, CoY = CoY,
                                             seed = i,
                                             max_itr = 60,
                                             max_tot.itr = 120,
                                             useY = TRUE))))

  betas <- mean(unlist(fit1$res_Beta$Beta))
  mus <- mean(unlist(fit1$res_Mu))
  sigma <- mean(unlist(fit1$res_Sigma))
  Gamma <- mean(numeric_leaf(fit1$res_Gamma))

  # check estimates are numerically stable
  expect_true(all(is.finite(numeric_leaf(fit1$res_Beta))))
  expect_true(all(is.finite(numeric_leaf(fit1$res_Mu))))
  expect_true(all(is.finite(numeric_leaf(fit1$res_Sigma))))
  expect_true(all(is.finite(numeric_leaf(fit1$res_Gamma))))
  expect_lt(abs(betas), 1)
  expect_lt(abs(mus), 1)
  expect_lt(abs(sigma), 1)
  expect_lt(abs(Gamma), 2)

  expect_equal(class(fit1), "lucid_serial")

  Z <- list(Z1 = Z1, list(Z2 = Z2, Z3 = Z3), Z4 = Z4, Z5 = Z5)
  suppressWarnings(invisible(capture.output(fit2 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2,list(2,2),2,2),
                                             CoG = CoG, CoY = CoY,
                                             lucid_model = "serial",
                                             family = "binary",
                                             init_omic.data.model = "VVV",
                                             seed = i,
                                             max_itr = 60,
                                             max_tot.itr = 120,
                                             useY = TRUE))))

  betas <- mean(unlist(fit2$res_Beta))
  mus <- mean(unlist(fit2$res_Mu))
  sigma <- mean(unlist(fit2$res_Sigma))
  Gamma <- mean(numeric_leaf(fit2$res_Gamma))

  # check estimates are numerically stable
  expect_true(all(is.finite(numeric_leaf(fit2$res_Beta))))
  expect_true(all(is.finite(numeric_leaf(fit2$res_Mu))))
  expect_true(all(is.finite(numeric_leaf(fit2$res_Sigma))))
  expect_true(all(is.finite(numeric_leaf(fit2$res_Gamma))))
  expect_lt(abs(betas), 1)
  expect_lt(abs(mus), 1)
  expect_lt(abs(sigma), 1)
  expect_lt(abs(Gamma), 2)

  expect_equal(class(fit2), "lucid_serial")


})

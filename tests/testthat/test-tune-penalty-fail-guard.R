skip_on_cran()

test_that("tune_lucid early stops cleanly when all grid fits fail", {
  set.seed(11)
  N <- 30
  G <- matrix(rnorm(N), nrow = N) # single exposure
  Z <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  expect_error(
    tune_lucid(
      G = G,
      Z = Z,
      Y = Y,
      lucid_model = "early",
      family = "normal",
      K = 2:3,
      Rho_G = c(0.1, 0.2), # requires >=2 exposure vars, so all fits fail
      Rho_Z_Mu = 0,
      Rho_Z_Cov = 0,
      max_itr = 5,
      max_tot.itr = 20,
      tol = 1e-1,
      seed = 11
    ),
    "fails to converge"
  )
})

test_that("tune_lucid parallel stops cleanly when all grid fits fail", {
  set.seed(12)
  N <- 30
  G <- matrix(rnorm(N), nrow = N) # single exposure
  Z1 <- matrix(rnorm(N * 3), nrow = N)
  Z2 <- matrix(rnorm(N * 3), nrow = N)
  Y <- rnorm(N)

  expect_error(
    tune_lucid(
      G = G,
      Z = list(Z1, Z2),
      Y = Y,
      lucid_model = "parallel",
      family = "normal",
      K = list(2:3, 2),
      Rho_G = c(0.1, 0.2), # requires >=2 exposure vars, so all fits fail
      Rho_Z_Mu = 0,
      Rho_Z_Cov = 0,
      max_itr = 5,
      max_tot.itr = 20,
      tol = 1e-1,
      seed = 12
    ),
    "fails to converge"
  )
})

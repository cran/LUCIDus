# Targeted regressions for internal workhorse utilities (P0/P1 fixes)

test_that("lucid_par_early forwards penalties and extracts exposure coefficients without intercept", {
  set.seed(123)
  n <- 12
  dat <- data.frame(
    g1 = rnorm(n),
    g2 = rnorm(n),
    z1 = rnorm(n),
    y = rnorm(n),
    cog1 = rnorm(n),
    coy1 = rnorm(n)
  )

  model <- list(
    K = 2,
    family = "normal",
    init_omic.data.model = "EEV",
    init_impute = "lod",
    init_par = "random",
    em_control = list(tol = 1e-3, max_itr = 5, max_tot.itr = 20),
    Rho = list(Rho_G = 0.12, Rho_Z_Mu = 0.34, Rho_Z_Cov = 0.56),
    res_Gamma = list(beta = c("(Intercept)" = 0, LC2 = 1))
  )

  fake_fit <- structure(
    list(
      res_Beta = matrix(
        c(
          0, 0, 0, 0,
          9, 10, 20, 30
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(NULL, c("intercept", "g1", "g2", "cog1"))
      ),
      res_Mu = matrix(c(0.1, -0.2), nrow = 2, dimnames = list(NULL, "z1")),
      res_Gamma = list(beta = c("(Intercept)" = 0.5, LC2 = 1.5))
    ),
    class = "early_lucid"
  )

  captured <- NULL
  testthat::local_mocked_bindings(
    est_lucid = function(...) {
      captured <<- list(...)
      fake_fit
    },
    .package = "LUCIDus"
  )

  pb <- progress::progress_bar$new(total = 2, show_after = 9999)
  pars <- suppressWarnings(
    lucid_par_early(
      data = dat,
      indices = seq_len(n),
      model = model,
      dimG = 2,
      dimZ = 1,
      dimCoY = 1,
      dimCoG = 1,
      prog = pb
    )
  )

  expect_equal(captured$Rho_G, model$Rho$Rho_G)
  expect_equal(captured$Rho_Z_Mu, model$Rho$Rho_Z_Mu)
  expect_equal(captured$Rho_Z_Cov, model$Rho$Rho_Z_Cov)

  expect_equal(unname(pars[1:2]), c(10, 20))
  expect_equal(names(pars)[1:2], c("g1.cluster2", "g2.cluster2"))
})

test_that("initialize_Delta binary with 5 layers uses layer 5 responsibilities", {
  set.seed(1001)
  n <- 120
  K <- rep(2, 5)

  mk_z <- function() {
    p2 <- plogis(rnorm(n))
    cbind(1 - p2, p2)
  }
  z <- lapply(seq_len(5), function(i) mk_z())

  # Outcome depends heavily on layer-5 second cluster responsibility
  y_prob <- plogis(-1 + 4 * (z[[5]][, 2] - 0.5))
  Y <- rbinom(n, size = 1, prob = y_prob)

  delta_a <- initialize_Delta(K = K, CoY = NULL, family = "binary", z = z, Y = Y)
  z_mod <- z
  z_mod[[5]] <- cbind(z[[5]][, 2], z[[5]][, 1]) # swap layer-5 responsibilities
  delta_b <- initialize_Delta(K = K, CoY = NULL, family = "binary", z = z_mod, Y = Y)

  expect_gt(max(abs(as.numeric(delta_a$mu) - as.numeric(delta_b$mu))), 1e-8)
})

test_that("Estep_early fallback keeps row-wise variation when mvn path fails", {
  N <- 4
  K <- 2
  G <- matrix(0, nrow = N, ncol = 1)
  Z <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3), nrow = N, byrow = TRUE)
  Y <- matrix(0, nrow = N, ncol = 1)
  beta <- matrix(0, nrow = K, ncol = 2) # intercept + one G
  mu <- matrix(c(0, 0, 1, 1), nrow = K, byrow = TRUE)
  sigma <- array(0, dim = c(2, 2, K))
  sigma[, , 1] <- matrix(c(1, NA, NA, 1), nrow = 2, byrow = TRUE) # trigger try-error
  sigma[, , 2] <- diag(2)
  family.list <- normal(K = K, dimCoY = 0)

  ll <- suppressWarnings(
    Estep_early(
      beta = beta,
      mu = mu,
      sigma = sigma,
      gamma = NULL,
      G = G,
      Z = Z,
      Y = Y,
      family.list = family.list,
      K = K,
      N = N,
      useY = FALSE,
      ind.na = rep(1, N),
      itr = 2,
      dimCoY = 0,
      CoY = NULL
    )
  )

  expect_equal(dim(ll), c(N, K))
  expect_gt(length(unique(round(ll[, 1], 8))), 1)
})

test_that("f_XtoZ returns informative error path for invalid mvn evaluation", {
  Z <- matrix(rnorm(8), nrow = 4, ncol = 2)
  Mu <- matrix(c(NA, 0, 1, 1), nrow = 2, byrow = TRUE)
  Sigma <- array(0, dim = c(2, 2, 2))
  Sigma[, , 1] <- diag(2)
  Sigma[, , 2] <- diag(2)

  expect_error(
    f_XtoZ(Z = Z, Mu_matrix = Mu, Sigma_matrix = Sigma),
    "Error in cluster"
  )
})

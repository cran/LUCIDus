## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(LUCIDus)

## ----simulate-associated-data-------------------------------------------------
make_serial_six_parallel_data <- function(n = 60, pG = 5, pZ = 3, seed = 5252) {
  set.seed(seed)

  # Exposure block
  G <- matrix(rnorm(n * pG), nrow = n, ncol = pG)
  colnames(G) <- paste0("G", seq_len(pG))

  # Covariates linked to exposures (induced association)
  CoG <- cbind(
    cov1 = G[, 1] + rnorm(n, sd = 0.25),
    cov2 = G[, 2] - 0.5 * G[, 3] + rnorm(n, sd = 0.25)
  )
  CoY <- CoG

  # Serial latent signals
  eta <- matrix(0, nrow = n, ncol = 6)
  x <- matrix(0, nrow = n, ncol = 6)

  eta[, 1] <- 0.9 * G[, 1] - 0.7 * G[, 2] + 0.4 * CoG[, 1] + rnorm(n, sd = 0.5)
  x[, 1] <- as.numeric(eta[, 1] > median(eta[, 1]))

  for (s in 2:6) {
    eta[, s] <- 0.7 * as.numeric(scale(eta[, s - 1])) +
      0.5 * G[, 1] - 0.4 * G[, 3] + 0.35 * x[, s - 1] + rnorm(n, sd = 0.6)
    x[, s] <- as.numeric(eta[, s] > median(eta[, s]))
  }

  # 6 stages; each stage has 2 parallel layers
  Z <- vector("list", 6)
  names(Z) <- paste0("stage", seq_len(6))

  for (s in seq_len(6)) {
    shared <- rnorm(n, sd = 0.35)

    layer1 <- cbind(
      1.2 * x[, s] + 0.5 * G[, 1] + shared + rnorm(n, sd = 0.45),
      0.9 * x[, s] - 0.3 * G[, 2] + shared + rnorm(n, sd = 0.45),
      0.7 * x[, s] + 0.4 * G[, 4] + shared + rnorm(n, sd = 0.45)
    )

    layer2 <- cbind(
      -1.0 * x[, s] + 0.45 * G[, 2] + shared + rnorm(n, sd = 0.45),
      -0.8 * x[, s] - 0.35 * G[, 1] + shared + rnorm(n, sd = 0.45),
      -0.6 * x[, s] + 0.30 * G[, 5] + shared + rnorm(n, sd = 0.45)
    )

    colnames(layer1) <- paste0("s", s, "_L1_f", seq_len(pZ))
    colnames(layer2) <- paste0("s", s, "_L2_f", seq_len(pZ))

    Z[[s]] <- list(layer1 = layer1, layer2 = layer2)
  }

  # Outcome associated with final latent stage + exposure/covariate terms
  Y <- 0.8 * x[, 6] + 0.45 * G[, 1] - 0.25 * G[, 2] + 0.35 * CoY[, 2] + rnorm(n, sd = 0.7)

  list(G = G, Z = Z, Y = as.numeric(Y), CoG = CoG, CoY = CoY)
}

d <- make_serial_six_parallel_data(seed = 5252)
str(d, max.level = 2)

## ----monotone-missingness-----------------------------------------------------
apply_monotone_missingness <- function(d, pZ = 3) {
  listwise_n_by_stage <- c(0, 4, 8, 12, 16, 20)

  for (s in seq_len(6)) {
    if (listwise_n_by_stage[s] > 0) {
      miss_rows <- seq_len(listwise_n_by_stage[s])
      d$Z[[s]]$layer1[miss_rows, ] <- NA
      d$Z[[s]]$layer2[miss_rows, ] <- NA
    }

    spor_rows <- seq_len(min(nrow(d$G), s + 1))
    d$Z[[s]]$layer1[spor_rows, 1] <- NA
    d$Z[[s]]$layer2[spor_rows, pZ] <- NA
  }

  d$listwise_n_by_stage <- listwise_n_by_stage
  d
}

d_miss <- apply_monotone_missingness(d)
d_miss$listwise_n_by_stage

## ----fit-serial-6parallel, warning=FALSE--------------------------------------
K6 <- replicate(6, list(2, 2), simplify = FALSE)

set.seed(5252)
serial_fit <- estimate_lucid(
  lucid_model = "serial",
  G = d_miss$G,
  Z = d_miss$Z,
  Y = d_miss$Y,
  CoG = d_miss$CoG,
  CoY = d_miss$CoY,
  family = "normal",
  K = K6,
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 6,
  max_tot.itr = 36,
  tol = 1e-2,
  seed = 5252,
  verbose = FALSE
)

class(serial_fit)
length(serial_fit$submodel)
sapply(serial_fit$submodel, class)

## ----fit-serial-verbose-demo, warning=FALSE-----------------------------------
K2 <- replicate(2, list(2, 2), simplify = FALSE)

# Smaller two-stage subset to keep verbose demo fast.
n_demo <- 30
d_demo <- list(
  G = d_miss$G[1:n_demo, , drop = FALSE],
  Z = lapply(d_miss$Z[1:2], function(stage) {
    list(
      layer1 = stage$layer1[1:n_demo, , drop = FALSE],
      layer2 = stage$layer2[1:n_demo, , drop = FALSE]
    )
  }),
  Y = d_miss$Y[1:n_demo],
  CoG = d_miss$CoG[1:n_demo, , drop = FALSE],
  CoY = d_miss$CoY[1:n_demo, , drop = FALSE]
)

set.seed(5251)
serial_fit_verbose <- estimate_lucid(
  lucid_model = "serial",
  G = d_demo$G,
  Z = d_demo$Z,
  Y = d_demo$Y,
  CoG = d_demo$CoG,
  CoY = d_demo$CoY,
  family = "normal",
  K = K2,
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 2,
  max_tot.itr = 8,
  tol = 1e-2,
  seed = 5251,
  verbose = TRUE
)

## ----check-missing-summary----------------------------------------------------
listwise_l1 <- vapply(serial_fit$missing_summary$stage, function(ms) {
  as.integer(ms$layer_summary$listwise_rows[1])
}, integer(1))

listwise_l2 <- vapply(serial_fit$missing_summary$stage, function(ms) {
  as.integer(ms$layer_summary$listwise_rows[2])
}, integer(1))

data.frame(
  stage = 1:6,
  expected_listwise = d_miss$listwise_n_by_stage,
  observed_layer1 = listwise_l1,
  observed_layer2 = listwise_l2
)

stopifnot(all(listwise_l1 == d_miss$listwise_n_by_stage))
stopifnot(all(listwise_l2 == d_miss$listwise_n_by_stage))
stopifnot(all(diff(listwise_l1) >= 0))

## ----check-final-stage--------------------------------------------------------
last_stage <- serial_fit$submodel[[6]]

all(is.finite(last_stage$inclusion.p[[1]]))
all(is.finite(last_stage$inclusion.p[[2]]))

stopifnot(all(is.finite(last_stage$inclusion.p[[1]])))
stopifnot(all(is.finite(last_stage$inclusion.p[[2]])))

## ----print-summary, warning=FALSE---------------------------------------------
summary(serial_fit)

## ----serial-bootstrap, warning=FALSE------------------------------------------
set.seed(5253)

serial_boot <- boot_lucid(
  G = d_miss$G,
  Z = d_miss$Z,
  Y = d_miss$Y,
  CoG = d_miss$CoG,
  CoY = d_miss$CoY,
  model = serial_fit,
  R = 2,
  conf = 0.90
)

# Stage-wise bootstrap objects are returned under $stage
length(serial_boot$stage)
names(serial_boot$stage[[1]])

## ----serial-summary-bootstrap, warning=FALSE----------------------------------
summary(serial_fit, boot.se = serial_boot)

## -----------------------------------------------------------------------------
sessionInfo()


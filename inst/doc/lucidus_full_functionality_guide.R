## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(LUCIDus)

## ----build-demo-data----------------------------------------------------------
make_demo_data <- function(n = 80, pG = 6, pZ = 4, seed = 20260309) {
  set.seed(seed)

  # Exposures
  G <- matrix(rnorm(n * pG), nrow = n, ncol = pG)
  colnames(G) <- paste0("G", seq_len(pG))

  # Covariates associated with exposures
  CoG <- cbind(
    cov_g1 = G[, 1] + 0.2 * G[, 2] + rnorm(n, sd = 0.2),
    cov_g2 = -0.3 * G[, 3] + 0.4 * G[, 4] + rnorm(n, sd = 0.25)
  )
  CoY <- CoG

  # Latent cluster driver
  lin <- 1.0 * G[, 1] - 0.8 * G[, 2] + 0.4 * CoG[, 1]
  prob_x <- plogis(lin)
  x <- rbinom(n, size = 1, prob = prob_x)

  # Layer 1 and 2 for parallel stage
  Z1 <- cbind(
    1.2 * x + 0.5 * G[, 1] + rnorm(n, sd = 0.5),
    1.0 * x - 0.4 * G[, 2] + rnorm(n, sd = 0.5),
    0.8 * x + 0.3 * G[, 3] + rnorm(n, sd = 0.5),
    0.6 * x + 0.2 * G[, 4] + rnorm(n, sd = 0.5)
  )

  Z2 <- cbind(
    -1.1 * x + 0.45 * G[, 2] + rnorm(n, sd = 0.5),
    -0.9 * x - 0.35 * G[, 1] + rnorm(n, sd = 0.5),
    -0.7 * x + 0.25 * G[, 5] + rnorm(n, sd = 0.5),
    -0.5 * x + 0.20 * G[, 6] + rnorm(n, sd = 0.5)
  )

  # Layer 3 for serial second stage (early stage)
  Z3 <- cbind(
    0.9 * x + 0.30 * G[, 1] + rnorm(n, sd = 0.55),
    0.7 * x - 0.25 * G[, 3] + rnorm(n, sd = 0.55),
    -0.8 * x + 0.20 * G[, 4] + rnorm(n, sd = 0.55),
    -0.6 * x + 0.15 * G[, 6] + rnorm(n, sd = 0.55)
  )

  colnames(Z1) <- paste0("Z1_f", seq_len(pZ))
  colnames(Z2) <- paste0("Z2_f", seq_len(pZ))
  colnames(Z3) <- paste0("Z3_f", seq_len(pZ))

  # Outcomes
  Y_normal <- 1.1 * x + 0.5 * G[, 1] - 0.25 * G[, 3] + 0.35 * CoY[, 2] + rnorm(n, sd = 0.7)
  Y_binary <- rbinom(n, size = 1, prob = plogis(-0.2 + 1.0 * x + 0.35 * G[, 1] - 0.2 * CoY[, 1]))

  # Structures for each model
  Z_parallel <- list(layer1 = Z1, layer2 = Z2)
  Z_early <- cbind(Z1, Z2)
  Z_serial_mixed <- list(list(layer1 = Z1, layer2 = Z2), Z3)

  list(
    G = G,
    CoG = CoG,
    CoY = CoY,
    Y_normal = as.numeric(Y_normal),
    Y_binary = as.numeric(Y_binary),
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3,
    Z_parallel = Z_parallel,
    Z_early = Z_early,
    Z_serial_mixed = Z_serial_mixed
  )
}

d <- make_demo_data()

## ----add-missingness----------------------------------------------------------
# Early matrix missingness
Z_early_miss <- d$Z_early
Z_early_miss[1, ] <- NA          # listwise
Z_early_miss[2:4, 1] <- NA       # sporadic block
Z_early_miss[5, 3] <- NA         # sporadic cell

# Parallel list missingness
Z_parallel_miss <- d$Z_parallel
Z_parallel_miss[[1]][1, ] <- NA  # listwise on layer1
Z_parallel_miss[[2]][2, 2] <- NA # sporadic on layer2

# Serial mixed missingness (stage1 parallel + stage2 early)
Z_serial_miss <- list(
  list(
    layer1 = Z_parallel_miss[[1]],
    layer2 = Z_parallel_miss[[2]]
  ),
  {tmp <- d$Z3; tmp[3, ] <- NA; tmp}
)

## ----missing-diagnostics------------------------------------------------------
miss_info_early <- analyze_missing_pattern(Z_early_miss)
na_early <- check_na(Z_early_miss, lucid_model = "early")
na_parallel <- check_na(Z_parallel_miss, lucid_model = "parallel")

miss_info_early$total_missing
table(na_early$indicator_na)
na_parallel$cross_layer_summary

## ----impute-helpers-----------------------------------------------------------
# Use a sub-matrix for concise display
orig_sub <- Z_early_miss[, 1:4, drop = FALSE]
imputed_sub <- safe_impute(orig_sub, method = "mean")
quality_sub <- check_imputation_quality(orig_sub, imputed_sub)

quality_sub$overall_quality

## ----fit-early, warning=FALSE-------------------------------------------------
set.seed(101)
fit_early <- estimate_lucid(
  lucid_model = "early",
  G = d$G,
  Z = Z_early_miss,
  Y = d$Y_normal,
  CoG = d$CoG,
  CoY = d$CoY,
  family = "normal",
  K = 2,
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 8,
  max_tot.itr = 30,
  tol = 1e-2,
  seed = 101,
  verbose = FALSE
)

class(fit_early)

## ----fit-early-structure------------------------------------------------------
names(fit_early)
cat("likelihood:", fit_early$likelihood, "\n")
cat("selected G:", sum(get_selected_G(fit_early)), "of",
    length(get_selected_G(fit_early)), "\n")
cat("selected Z:", sum(get_selected_Z(fit_early)), "of",
    length(get_selected_Z(fit_early)), "\n")

## ----fit-parallel, warning=FALSE----------------------------------------------
set.seed(102)
fit_parallel <- estimate_lucid(
  lucid_model = "parallel",
  G = d$G,
  Z = Z_parallel_miss,
  Y = d$Y_normal,
  CoG = d$CoG,
  CoY = d$CoY,
  family = "normal",
  K = c(2, 2),
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 8,
  max_tot.itr = 30,
  tol = 1e-2,
  seed = 102,
  verbose = FALSE
)

class(fit_parallel)

## ----fit-parallel-structure---------------------------------------------------
names(fit_parallel)
cat("N:", fit_parallel$N, "\n")
cat("selected G per layer:", paste(sapply(seq_along(fit_parallel$K), function(i)
                                     sum(get_selected_G(fit_parallel, layer = i))),
                                   collapse = ", "), "\n")

sel_z_parallel <- get_selected_Z(fit_parallel)
cat("selected Z per layer:", paste(sapply(sel_z_parallel, sum), collapse = ", "), "\n")

## ----fit-serial, warning=FALSE------------------------------------------------
set.seed(103)
fit_serial <- estimate_lucid(
  lucid_model = "serial",
  G = d$G,
  Z = Z_serial_miss,
  Y = d$Y_normal,
  CoG = d$CoG,
  CoY = d$CoY,
  family = "normal",
  K = list(list(2, 2), 2),
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 8,
  max_tot.itr = 36,
  tol = 1e-2,
  seed = 103,
  verbose = FALSE
)

class(fit_serial)
length(fit_serial$submodel)

## ----fit-serial-structure-----------------------------------------------------
names(fit_serial)
cat("top-level likelihood (sum over stages):", fit_serial$likelihood, "\n")
cat("top-level select is stage 1's select:",
    identical(fit_serial$select, fit_serial$submodel[[1]]$select), "\n")

sel_z_serial <- get_selected_Z(fit_serial)
cat("stage 1 (parallel) selected Z per layer:",
    paste(sapply(sel_z_serial[[1]], sum), collapse = ", "), "\n")
cat("stage 2 (early) selected Z:", sum(sel_z_serial[[2]]), "\n")

## ----fit-verbose-demos, warning=FALSE-----------------------------------------
n_demo <- 30
G_demo <- d$G[1:n_demo, , drop = FALSE]
Y_demo <- d$Y_normal[1:n_demo]
CoG_demo <- d$CoG[1:n_demo, , drop = FALSE]
CoY_demo <- d$CoY[1:n_demo, , drop = FALSE]
Z_early_demo <- d$Z_early[1:n_demo, , drop = FALSE]
Z_parallel_demo <- lapply(d$Z_parallel, function(z) z[1:n_demo, , drop = FALSE])
Z_serial_demo <- list(
  lapply(d$Z_parallel, function(z) z[1:n_demo, , drop = FALSE]),
  d$Z3[1:n_demo, , drop = FALSE]
)

set.seed(111)
fit_early_verbose <- estimate_lucid(
  lucid_model = "early",
  G = G_demo,
  Z = Z_early_demo,
  Y = Y_demo,
  CoG = CoG_demo,
  CoY = CoY_demo,
  family = "normal",
  K = 2,
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 2,
  max_tot.itr = 8,
  tol = 1e-2,
  seed = 111,
  verbose = TRUE
)

set.seed(112)
fit_parallel_verbose <- estimate_lucid(
  lucid_model = "parallel",
  G = G_demo,
  Z = Z_parallel_demo,
  Y = Y_demo,
  CoG = CoG_demo,
  CoY = CoY_demo,
  family = "normal",
  K = c(2, 2),
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 2,
  max_tot.itr = 8,
  tol = 1e-2,
  seed = 112,
  verbose = TRUE
)

set.seed(113)
fit_serial_verbose <- estimate_lucid(
  lucid_model = "serial",
  G = G_demo,
  Z = Z_serial_demo,
  Y = Y_demo,
  CoG = CoG_demo,
  CoY = CoY_demo,
  family = "normal",
  K = list(list(2, 2), 2),
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 2,
  max_tot.itr = 10,
  tol = 1e-2,
  seed = 113,
  verbose = TRUE
)

## ----tune-lucid, warning=FALSE------------------------------------------------
set.seed(104)

tune_early <- tune_lucid(
  G = d$G,
  Z = d$Z_early,
  Y = d$Y_normal,
  CoG = d$CoG,
  CoY = d$CoY,
  family = "normal",
  lucid_model = "early",
  K = 2:3,
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 6,
  max_tot.itr = 24,
  seed = 104
)

head(tune_early$tune_list)

## ----lucid-wrapper, warning=FALSE---------------------------------------------
set.seed(105)

fit_lucid_wrapper <- lucid(
  G = d$G,
  Z = d$Z_early,
  Y = d$Y_normal,
  CoG = d$CoG,
  CoY = d$CoY,
  family = "normal",
  lucid_model = "early",
  K = 2:3,
  Rho_G = 0,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 6,
  max_tot.itr = 24,
  seed = 105
)

class(fit_lucid_wrapper)

## ----summaries, warning=FALSE-------------------------------------------------
summary(fit_early)
summary(fit_parallel)
summary(fit_serial)

## ----predict-standard, warning=FALSE------------------------------------------
# Use lightweight no-covariate fits for robust prediction demo.
set.seed(205)
fit_early_pred <- estimate_lucid(
  lucid_model = "early",
  G = d$G,
  Z = d$Z_early,
  Y = d$Y_normal,
  family = "normal",
  K = 2,
  max_itr = 6,
  max_tot.itr = 20,
  tol = 1e-2,
  seed = 205
)

set.seed(206)
fit_parallel_pred <- estimate_lucid(
  lucid_model = "parallel",
  G = d$G,
  Z = d$Z_parallel,
  Y = d$Y_normal,
  family = "normal",
  K = c(2, 2),
  max_itr = 6,
  max_tot.itr = 20,
  tol = 1e-2,
  seed = 206
)

set.seed(207)
fit_serial_pred <- estimate_lucid(
  lucid_model = "serial",
  G = d$G,
  Z = d$Z_serial_mixed,
  Y = d$Y_normal,
  family = "normal",
  K = list(list(2, 2), 2),
  max_itr = 6,
  max_tot.itr = 24,
  tol = 1e-2,
  seed = 207
)

pred_early <- predict_lucid(
  model = fit_early_pred,
  G = d$G,
  Z = d$Z_early,
  Y = d$Y_normal
)

pred_parallel <- predict_lucid(
  model = fit_parallel_pred,
  G = d$G,
  Z = d$Z_parallel,
  Y = d$Y_normal
)

pred_serial <- predict_lucid(
  model = fit_serial_pred,
  G = d$G,
  Z = d$Z_serial_mixed,
  Y = d$Y_normal
)

# Cluster assignment. pred.x is a vector for early and a list -- by layer for
# parallel, by stage for serial -- so the blocks are summarised in one place.
# pred.x nests: a vector for early, a list by layer for parallel, and for
# serial a list by stage whose elements are themselves lists when that stage is
# a parallel submodel. Flatten to the leaves so each cluster variable is counted
# on its own -- pooling a parallel stage's layers would report twice as many
# assignments as there are subjects.
flatten_blocks <- function(x, path = "") {
  if (!is.list(x)) return(stats::setNames(list(as.numeric(x)), path))
  out <- list()
  for (i in seq_along(x)) {
    nm <- if (nzchar(path)) paste0(path, ".", i) else as.character(i)
    out <- c(out, flatten_blocks(x[[i]], nm))
  }
  out
}

cluster_sizes <- function(pred_x, label) {
  blocks <- flatten_blocks(pred_x)
  nms <- names(blocks)
  # index by position: the single block of an early model is named "", and
  # blocks[[""]] does not select anything.
  do.call(rbind, lapply(seq_along(blocks), function(i) {
    tb <- table(factor(blocks[[i]]))
    data.frame(model = label,
               block = if (!nzchar(nms[i])) "-" else nms[i],
               cluster = names(tb), n = as.integer(tb), row.names = NULL)
  }))
}

rbind(
  cluster_sizes(pred_early$pred.x,    "early"),
  cluster_sizes(pred_parallel$pred.x, "parallel"),
  cluster_sizes(pred_serial$pred.x,   "serial")
)

## ----cluster-assignment-extractor---------------------------------------------
identical(as.numeric(get_cluster_assignment(fit_early_pred)),
          as.numeric(pred_early$pred.x))

## ----predict-standard-y-------------------------------------------------------
outcome_summary <- function(pred_y, label) {
  v <- as.numeric(unlist(pred_y))
  data.frame(
    model  = label,
    n      = length(v),
    mean   = round(mean(v), 3),
    sd     = round(stats::sd(v), 3),
    min    = round(min(v), 3),
    median = round(stats::median(v), 3),
    max    = round(max(v), 3),
    row.names = NULL
  )
}

rbind(
  outcome_summary(pred_early$pred.y,    "early"),
  outcome_summary(pred_parallel$pred.y, "parallel"),
  outcome_summary(pred_serial$pred.y,   "serial")
)

## ----predict-standard-vs-observed---------------------------------------------
rbind(
  outcome_summary(pred_early$pred.y, "predicted (early)"),
  outcome_summary(d$Y_normal,        "observed")
)

## ----predict-gcomp, warning=FALSE---------------------------------------------
pred_early_g <- predict_lucid(
  model = fit_early_pred,
  G = d$G,
  Z = NULL,
  Y = NULL,
  g_computation = TRUE
)

pred_parallel_g <- try(
  predict_lucid(
    model = fit_parallel_pred,
    G = d$G,
    Z = NULL,
    Y = NULL,
    g_computation = TRUE
  ),
  silent = TRUE
)

names(pred_early_g)
if (inherits(pred_parallel_g, "try-error")) {
  "parallel g_computation returned a try-error on this demo object; code pattern is shown above."
} else {
  names(pred_parallel_g)
}

pred_serial_g <- predict_lucid(
  model = fit_serial_pred,
  G = d$G,
  Z = NULL,
  Y = NULL,
  g_computation = TRUE
)

names(pred_serial_g)
length(pred_serial_g$pred.z)

## ----predict-modes------------------------------------------------------------
pred_unsup <- predict_lucid(
  model = fit_early_pred,
  G = d$G,
  Z = d$Z_early
)

missing_z <- try(
  predict_lucid(
    model = fit_early_pred,
    G = d$G,
    Z = NULL,
    Y = d$Y_normal
  ),
  silent = TRUE
)

data.frame(
  mode = c("Y omitted (unsupervised)", "Z omitted, no g-computation"),
  result = c(
    paste(names(pred_unsup), collapse = ", "),
    if (inherits(missing_z, "try-error")) "error, as documented" else "unexpectedly succeeded"
  )
)

cat(if (inherits(missing_z, "try-error")) attr(missing_z, "condition")$message else "")

## ----bootstrap-all, warning=FALSE---------------------------------------------
set.seed(106)

boot_early <- boot_lucid(
  G = d$G,
  Z = Z_early_miss,
  Y = d$Y_normal,
  CoG = d$CoG,
  CoY = d$CoY,
  model = fit_early,
  R = 2,
  conf = 0.9
)

boot_parallel <- boot_lucid(
  G = d$G,
  Z = Z_parallel_miss,
  Y = d$Y_normal,
  CoG = d$CoG,
  CoY = d$CoY,
  model = fit_parallel,
  R = 2,
  conf = 0.9
)

boot_serial <- boot_lucid(
  G = d$G,
  Z = Z_serial_miss,
  Y = d$Y_normal,
  CoG = d$CoG,
  CoY = d$CoY,
  model = fit_serial,
  R = 2,
  conf = 0.9
)

summary(fit_early, boot.se = boot_early)
summary(fit_parallel, boot.se = boot_parallel)
summary(fit_serial, boot.se = boot_serial)

## ----plotting, warning=FALSE--------------------------------------------------
plot_early <- plot(fit_early)
plot_parallel <- try(plot(fit_parallel), silent = TRUE)
plot_serial <- try(plot(fit_serial), silent = TRUE)

class(plot_early)
if (inherits(plot_parallel, "try-error")) {
  "plot(fit_parallel) returned try-error in this build (parallel plot is under development)."
} else {
  class(plot_parallel)
}
if (inherits(plot_serial, "try-error")) {
  "plot(fit_serial) returned try-error in this build (serial plot is under development)."
} else {
  class(plot_serial)
}

## ----profile-heatmap-early, fig.width = 6, fig.height = 4.5-------------------
prof_early <- plot_cluster_omic_profile(fit_early, top_n = 8)
names(prof_early)
prof_early[[1]]

## ----profile-bar-early, fig.width = 6.5, fig.height = 4.5---------------------
plot_cluster_omic_profile(fit_early, type = "bar", top_n = 8)[[1]]

## ----profile-parallel, fig.width = 6.5, fig.height = 4.5----------------------
prof_parallel <- plot_cluster_omic_profile(fit_parallel, top_n = 8)
names(prof_parallel)
prof_parallel[[1]]

## ----profile-parallel-bar, fig.width = 6.5, fig.height = 4.5------------------
plot_cluster_omic_profile(fit_parallel, type = "bar", top_n = 8)[[2]]

## ----profile-serial, fig.width = 6.5, fig.height = 4.5------------------------
prof_serial <- plot_cluster_omic_profile(fit_serial, top_n = 8)
names(prof_serial)
prof_serial[[1]]

## ----profile-serial-bar, fig.width = 6.5, fig.height = 4.5--------------------
plot_cluster_omic_profile(fit_serial, type = "bar", top_n = 8)[[2]]

## ----top-omics-features-------------------------------------------------------
top_early <- get_top_omics_features(fit_early, top_n = 5)
top_early

## ----binary-example, warning=FALSE--------------------------------------------
set.seed(107)
fit_early_binary <- estimate_lucid(
  lucid_model = "early",
  G = d$G,
  Z = d$Z_early,
  Y = d$Y_binary,
  CoG = d$CoG,
  CoY = d$CoY,
  family = "binary",
  K = 2,
  max_itr = 8,
  max_tot.itr = 30,
  tol = 1e-2,
  seed = 107
)

pred_binary_prob <- predict_lucid(
  model = fit_early_binary,
  G = d$G,
  Z = d$Z_early,
  Y = d$Y_binary,
  CoG = d$CoG,
  CoY = d$CoY,
  response = FALSE
)

range(pred_binary_prob$pred.y)

## -----------------------------------------------------------------------------
sessionInfo()


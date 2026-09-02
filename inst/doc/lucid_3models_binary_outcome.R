## ----setup, message=FALSE, warning=FALSE--------------------------------------
# Keep knitting on error, so that one failing step reports itself and the rest of
# the tutorial still runs. The status table in section 13 records what happened.
knitr::opts_chunk$set(error = TRUE)

# Lightweight registry so the document can verify itself rather than relying on
# the reader to notice a missing output.
.reg <- new.env(parent = emptyenv()); .reg$rows <- list()
check_obj <- function(name, expected_class = NULL, section = "") {
  ok <- exists(name, envir = globalenv())
  cls <- if (ok) class(get(name, envir = globalenv()))[1] else NA_character_
  status <- if (!ok) "MISSING"
            else if (!is.null(expected_class) && !identical(cls, expected_class)) "unexpected class"
            else "ok"
  .reg$rows[[length(.reg$rows) + 1L]] <-
    data.frame(section = section, object = name, class = cls,
               status = status, stringsAsFactors = FALSE)
  invisible(NULL)
}

library(LUCIDus)

# The HELIX simulation object bundled with the package.
data(simulated_HELIX_data)

## ----data-prep----------------------------------------------------------------
# Use a smaller subset for vignette speed while preserving model behavior.
idx <- 1:90
ph <- simulated_HELIX_data$phenotype[idx, ]
n <- nrow(ph)

set.seed(2026)

# ---------------------------------------------------------------------------
# A tutorial dataset with a KNOWN answer.
#
# The HELIX omics matrices are real simulated data with their own structure, and
# the exposures shipped with them have no relationship to it. That is fine for
# demonstrating that code runs, but it makes feature selection impossible to
# judge: there is no right answer to compare against. So we plant one.
#
# The generating story, which is the DAG LUCID assumes:
#
#     causal exposures  ->  latent subgroup  ->  omics profile
#                                            ->  outcome
#
# Three exposures carry the subgroup signal with graded strength; six are pure
# noise. Half the features of each omics layer are shifted by subgroup
# membership; the rest are left as they came. Selection therefore has an
# unambiguous target, and the tutorial can check its answer instead of asserting
# it.
# ---------------------------------------------------------------------------

# The true latent subgroup. Retained so every selection claim below can be
# checked against it.
x_true <- rbinom(n, 1, 0.5)

# Exposures. g_causal_* predict subgroup membership; g_noise_* do not.
# Effect sizes are deliberately moderate. Stronger exposures make selection
# look better but drive the G -> X model to saturation, where every subject sits
# at posterior probability 1 and no counterfactual shift can move anything --
# which would make the g-computation demonstration in the continuous-outcome
# companion vignette vacuous.
G <- cbind(
  g_causal_1 =  1.0 * (x_true - 0.5) + rnorm(n, sd = 0.8),   # strongest
  g_causal_2 = -0.8 * (x_true - 0.5) + rnorm(n, sd = 0.8),   # moderate, negative
  g_causal_3 =  0.6 * (x_true - 0.5) + rnorm(n, sd = 0.8),   # weakest
  g_noise_1 = rnorm(n), g_noise_2 = rnorm(n), g_noise_3 = rnorm(n),
  g_noise_4 = rnorm(n), g_noise_5 = rnorm(n), g_noise_6 = rnorm(n)
)
G <- as.matrix(scale(G))

causal_exposures <- c("g_causal_1", "g_causal_2", "g_causal_3")

# Exposure penalty used throughout. The continuous-outcome companion vignette
# shows what this value recovers and sweeps it, along with the omics penalty,
# separately.
RHO_G <- 0.05

# Covariates for G->X (CoG) and X->Y (CoY).
# Here we use age-related and sex covariates from phenotype.
CoG <- cbind(
  hs_child_age_yrs_None = as.numeric(ph$hs_child_age_yrs_None),
  sex_male = as.numeric(ph$e3_sex_None == "male")
)
CoY <- CoG

# Two outcomes on the SAME subjects, so the normal and binary results below are
# directly comparable: the only thing that changes between them is the outcome
# model, not the sample, the omics, or the injected missingness.
#
# Continuous outcome: the real CK-18 measurement, plus a subgroup effect so the
# cluster -> outcome arm of the model has something to estimate.
Y <- as.numeric(ph$ck18_scaled) + 1.2 * x_true

# Binary outcome: median split. The median is used rather than a higher
# threshold because it splits these 90 subjects 45/45, and a balanced outcome
# gives the K = 2 outcome model the most to work with at this sample size.
Y_binary <- as.integer(Y > median(Y))
cat("binary outcome balance:\n"); print(table(Y_binary))

# Construct three omics layers and standardize each, then plant the subgroup
# signal in the first three features of every layer. The remaining seven per
# layer are left as they came and act as omics noise.
meth <- scale(simulated_HELIX_data$methylome[idx, 1:10, drop = FALSE])
tran <- scale(simulated_HELIX_data$transcriptome[idx, 1:10, drop = FALSE])
mir  <- scale(simulated_HELIX_data$miRNA[idx, 1:10, drop = FALSE])

signal_features <- 1:3
omics_shift <- 3.0
meth[, signal_features] <- meth[, signal_features] + omics_shift * x_true
tran[, signal_features] <- tran[, signal_features] - omics_shift * x_true
mir[,  signal_features] <- mir[,  signal_features] + omics_shift * x_true

# Column positions of the signal features once the layers are stacked for the
# early model, so selection can be scored against them later.
signal_cols_early <- c(signal_features,
                       ncol(meth) + signal_features,
                       ncol(meth) + ncol(tran) + signal_features)

# Early model uses one combined Z matrix.
Z_early <- cbind(meth, tran, mir)

# Parallel model uses list-of-layers.
Z_parallel <- list(methylome = meth, transcriptome = tran, miRNA = mir)

# Inject listwise + sporadic missingness for demonstration.
Z_early_miss <- Z_early
Z_early_miss[1, ] <- NA      # listwise row
Z_early_miss[2:4, 1] <- NA   # sporadic block
Z_early_miss[5, 3] <- NA     # sporadic cell

Z_parallel_miss <- Z_parallel
Z_parallel_miss[[1]][1, ] <- NA  # listwise in layer 1
Z_parallel_miss[[2]][2, 2] <- NA # sporadic in layer 2
Z_parallel_miss[[3]][3, 1] <- NA # sporadic in layer 3

# Quick structural sanity check.
str(list(
  G = G,
  CoG = CoG,
  CoY = CoY,
  Y = Y,
  Z_early = Z_early_miss,
  Z_parallel = Z_parallel_miss
), max.level = 1)

## ----helper-functions---------------------------------------------------------
# get_selected_G()/get_selected_Z() (from the package itself) already return a
# well-shaped, aligned logical mask straight from the fitted object -- no
# length mismatch is possible, since they derive it from the model's own
# recorded fields. The one thing left for a tutorial to decide is what to do
# if a penalty happened to deselect EVERY feature: refitting on zero columns
# would fail, so this keeps everything instead in that one edge case.
keep_or_all <- function(mask) if (any(mask, na.rm = TRUE)) mask else rep(TRUE, length(mask))

# Build selected-only inputs for early model.
prepare_early_selected_inputs <- function(fit_pen, G, Z) {
  list(
    G = as.matrix(G[, keep_or_all(get_selected_G(fit_pen)), drop = FALSE]),
    Z = as.matrix(Z[, keep_or_all(get_selected_Z(fit_pen)), drop = FALSE])
  )
}

# Build selected-only inputs for parallel model.
prepare_parallel_selected_inputs <- function(fit_pen, G, Z) {
  keep_g <- keep_or_all(get_selected_G(fit_pen))
  Z_sel <- lapply(seq_along(Z), function(i) {
    zi <- as.matrix(Z[[i]])
    zi[, keep_or_all(get_selected_Z(fit_pen, layer = i)), drop = FALSE]
  })
  names(Z_sel) <- names(Z)
  list(
    G = as.matrix(G[, keep_g, drop = FALSE]),
    Z = Z_sel
  )
}

# Serial stage>1 uses latent-cluster-derived "G" internally.
# We therefore subset stage-1 original G and each stage's Z where applicable.
prepare_serial_selected_inputs <- function(fit_pen, G, Z) {
  G_refit <- as.matrix(G)
  keep_g1 <- get_selected_G(fit_pen)
  if (length(keep_g1) == ncol(G_refit)) {
    G_refit <- G_refit[, keep_or_all(keep_g1), drop = FALSE]
  }

  selected_z <- get_selected_Z(fit_pen)
  Z_refit <- Z
  for (i in seq_along(fit_pen$submodel)) {
    sm <- fit_pen$submodel[[i]]
    if (inherits(sm, "early_lucid")) {
      zi <- as.matrix(Z_refit[[i]])
      Z_refit[[i]] <- zi[, keep_or_all(selected_z[[i]]), drop = FALSE]
    } else if (inherits(sm, "lucid_parallel")) {
      zi_list <- Z_refit[[i]]
      for (j in seq_along(zi_list)) {
        zij <- as.matrix(zi_list[[j]])
        zi_list[[j]] <- zij[, keep_or_all(selected_z[[i]][[j]]), drop = FALSE]
      }
      Z_refit[[i]] <- zi_list
    }
  }

  list(G = G_refit, Z = Z_refit)
}

# Zero-penalty refit, for any model type.
#
# The three model types previously had three byte-identical wrappers differing
# only in `lucid_model` and whether `useY` was forwarded; they are one function
# here. Everything about the model -- family, K, initialization, EM controls --
# is carried over from the screening fit, so the ONLY difference between the
# screening fit and this one is that the penalties are zero. That is what makes
# the refit estimates unshrunk and therefore suitable for bootstrap inference.
refit_selected <- function(model_type, fit_pen, inputs, Y,
                           CoG = NULL, CoY = NULL, seed = 1, verbose = FALSE) {
  args <- list(
    lucid_model = model_type,
    G = inputs$G,
    Z = inputs$Z,
    Y = Y,
    CoG = CoG,
    CoY = CoY,
    family = fit_pen$family,
    K = fit_pen$K,
    init_omic.data.model = fit_pen$init_omic.data.model,
    init_impute = fit_pen$init_impute,
    init_par = fit_pen$init_par,
    Rho_G = 0,
    Rho_Z_Mu = 0,
    Rho_Z_Cov = 0,
    max_itr = fit_pen$em_control$max_itr,
    max_tot.itr = fit_pen$em_control$max_tot.itr,
    tol = fit_pen$em_control$tol,
    seed = seed,
    verbose = verbose
  )
  # Every fitted class records useY, so it is carried over for all three model
  # types. The original three wrappers omitted it on the early path, which meant
  # an unsupervised screening fit would have been silently refitted supervised.
  args$useY <- fit_pen$useY
  do.call(estimate_lucid, args)
}

# Dispatcher for the three input-preparation helpers above.
prepare_selected_inputs <- function(model_type, fit_pen, G, Z) {
  switch(model_type,
    early    = prepare_early_selected_inputs(fit_pen, G, Z),
    parallel = prepare_parallel_selected_inputs(fit_pen, G, Z),
    serial   = prepare_serial_selected_inputs(fit_pen, G, Z),
    stop("unknown model_type: ", model_type)
  )
}


# Compact stage-wise feature-selection report for serial fits, built entirely
# from get_selected_G()/get_selected_Z() -- no per-stage dispatch of its own.
serial_selection_report <- function(fit_serial_pen) {
  selected_z <- get_selected_Z(fit_serial_pen)
  out <- vector("list", length(fit_serial_pen$submodel))
  for (i in seq_along(fit_serial_pen$submodel)) {
    sm <- fit_serial_pen$submodel[[i]]
    if (inherits(sm, "early_lucid")) {
      out[[i]] <- list(
        stage = i,
        model = "early",
        selected_G = if (i == 1) sum(get_selected_G(fit_serial_pen)) else NA,
        total_G = if (i == 1) length(get_selected_G(fit_serial_pen)) else NA,
        selected_Z = sum(selected_z[[i]]),
        total_Z = length(selected_z[[i]])
      )
    } else {
      out[[i]] <- list(
        stage = i,
        model = "parallel",
        selected_G = if (i == 1) sum(get_selected_G(fit_serial_pen)) else NA,
        total_G = if (i == 1) length(get_selected_G(fit_serial_pen)) else NA,
        selected_Z_by_layer = sapply(selected_z[[i]], sum),
        total_Z_by_layer = sapply(selected_z[[i]], length)
      )
    }
  }
  out
}

## ----early-binary-penalized, warning=FALSE------------------------------------
set.seed(1105)

early_pen_bin <- estimate_lucid(
  lucid_model = "early",
  G = G,
  Z = Z_early_miss,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  family = "binary",
  K = 2,
  Rho_G = RHO_G,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 15,
  max_tot.itr = 40,
  tol = 1e-2,
  seed = 1105,
  verbose = FALSE
)

summary(early_pen_bin)

## ----early-binary-refit, warning=FALSE----------------------------------------
set.seed(1106)

early_inputs_bin <- prepare_early_selected_inputs(early_pen_bin, G, Z_early_miss)

early_bin <- list(fit_pen = early_pen_bin, inputs = early_inputs_bin)
early_bin$fit_refit <- refit_selected(
  "early",
  fit_pen = early_pen_bin,
  inputs = early_inputs_bin,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  seed = 1106
)

summary(early_bin$fit_refit)

## ----early-binary-bootstrap, warning=FALSE------------------------------------
set.seed(1107)

early_bin$boot <- boot_lucid(
  G = early_inputs_bin$G,
  Z = early_inputs_bin$Z,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  model = early_bin$fit_refit,
  R = 30,
  conf = 0.90
)

summary(early_bin$fit_refit, boot.se = early_bin$boot)

## ----parallel-binary-penalized, warning=FALSE---------------------------------
set.seed(1205)

parallel_pen_bin <- estimate_lucid(
  lucid_model = "parallel",
  G = G,
  Z = Z_parallel_miss,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  family = "binary",
  K = c(2, 2, 2),
  Rho_G = RHO_G,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 15,
  max_tot.itr = 40,
  tol = 1e-2,
  seed = 1205,
  verbose = FALSE
)

summary(parallel_pen_bin)

## ----parallel-binary-refit, warning=FALSE-------------------------------------
set.seed(1206)

parallel_inputs_bin <- prepare_parallel_selected_inputs(parallel_pen_bin, G, Z_parallel_miss)

parallel_bin <- list(fit_pen = parallel_pen_bin, inputs = parallel_inputs_bin)
parallel_bin$fit_refit <- refit_selected(
  "parallel",
  fit_pen = parallel_pen_bin,
  inputs = parallel_inputs_bin,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  seed = 1206
)

summary(parallel_bin$fit_refit)

## ----parallel-binary-bootstrap, warning=FALSE---------------------------------
set.seed(1207)

parallel_bin$boot <- boot_lucid(
  G = parallel_inputs_bin$G,
  Z = parallel_inputs_bin$Z,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  model = parallel_bin$fit_refit,
  R = 30,
  conf = 0.90
)

summary(parallel_bin$fit_refit, boot.se = parallel_bin$boot)

## ----serial-all-early-binary-penalized, warning=FALSE-------------------------
# Serial structure: list of early-stage matrices.
Z_serial_all_early <- list(
  methylome = Z_parallel_miss[[1]],
  transcriptome = Z_parallel_miss[[2]],
  miRNA = Z_parallel_miss[[3]]
)

set.seed(1305)

serial_ae_pen_bin <- estimate_lucid(
  lucid_model = "serial",
  G = G,
  Z = Z_serial_all_early,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  family = "binary",
  K = list(2, 2, 2),
  Rho_G = RHO_G,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 15,
  max_tot.itr = 40,
  tol = 1e-2,
  seed = 1305,
  verbose = FALSE
)

summary(serial_ae_pen_bin)

## ----serial-all-early-binary-refit, warning=FALSE-----------------------------
set.seed(1306)

serial_ae_inputs_bin <- prepare_serial_selected_inputs(serial_ae_pen_bin, G, Z_serial_all_early)

serial_ae_bin <- list(fit_pen = serial_ae_pen_bin, inputs = serial_ae_inputs_bin)
serial_ae_bin$fit_refit <- refit_selected(
  "serial",
  fit_pen = serial_ae_pen_bin,
  inputs = serial_ae_inputs_bin,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  seed = 1306
)

summary(serial_ae_bin$fit_refit)

## ----serial-all-early-binary-bootstrap, warning=FALSE-------------------------
set.seed(1307)

serial_ae_bin$boot <- boot_lucid(
  G = serial_ae_inputs_bin$G,
  Z = serial_ae_inputs_bin$Z,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  model = serial_ae_bin$fit_refit,
  R = 30,
  conf = 0.90
)

summary(serial_ae_bin$fit_refit, boot.se = serial_ae_bin$boot)

## ----serial-mixed-binary-penalized, warning=FALSE-----------------------------
# Nested list signals a parallel submodel at stage 1, followed by early stage 2.
Z_serial_mixed <- list(
  list(
    methylome = Z_parallel_miss[[1]],
    transcriptome = Z_parallel_miss[[2]]
  ),
  miRNA = Z_parallel_miss[[3]]
)

set.seed(1405)

serial_mixed_pen_bin <- estimate_lucid(
  lucid_model = "serial",
  G = G,
  Z = Z_serial_mixed,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  family = "binary",
  K = list(list(2, 2), 2),
  Rho_G = RHO_G,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 15,
  max_tot.itr = 40,
  tol = 1e-2,
  seed = 1405,
  verbose = FALSE
)

summary(serial_mixed_pen_bin)

## ----serial-mixed-binary-refit, warning=FALSE---------------------------------
set.seed(1406)

serial_mixed_inputs_bin <- prepare_serial_selected_inputs(serial_mixed_pen_bin, G, Z_serial_mixed)

serial_mixed_bin <- list(fit_pen = serial_mixed_pen_bin, inputs = serial_mixed_inputs_bin)
serial_mixed_bin$fit_refit <- refit_selected(
  "serial",
  fit_pen = serial_mixed_pen_bin,
  inputs = serial_mixed_inputs_bin,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  seed = 1406
)

summary(serial_mixed_bin$fit_refit)

## ----serial-mixed-binary-bootstrap, warning=FALSE-----------------------------
set.seed(1407)

serial_mixed_bin$boot <- boot_lucid(
  G = serial_mixed_inputs_bin$G,
  Z = serial_mixed_inputs_bin$Z,
  Y = Y_binary,
  CoG = CoG,
  CoY = CoY,
  model = serial_mixed_bin$fit_refit,
  R = 30,
  conf = 0.90
)

summary(serial_mixed_bin$fit_refit, boot.se = serial_mixed_bin$boot)

## ----sankey-binary------------------------------------------------------------
sankey_early_bin <- plot(early_bin$fit_refit)
sankey_early_bin

## ----profile-early-binary, fig.width = 6.5, fig.height = 4.5, warning=FALSE----
prof_early_bin <- plot_cluster_omic_profile(early_bin$fit_refit, top_n = 10)
prof_early_bin[[1]]

## ----profile-parallel-binary, fig.width = 6.5, fig.height = 4.5, warning=FALSE----
prof_par_bin <- plot_cluster_omic_profile(
  parallel_bin$fit_refit,
  layer_names = c("methylome", "transcriptome", "miRNA"),
  top_n = 8
)
for (nm in names(prof_par_bin)) print(prof_par_bin[[nm]])

## ----profile-serial-binary, fig.width = 6.5, fig.height = 4.5, warning=FALSE----
prof_ser_bin <- plot_cluster_omic_profile(serial_ae_bin$fit_refit, top_n = 8)
for (nm in names(prof_ser_bin)) print(prof_ser_bin[[nm]])

## ----predict-binary-----------------------------------------------------------
pred_lab <- predict_lucid(model = early_bin$fit_refit,
                          G = early_bin$inputs$G, Z = early_bin$inputs$Z,
                          CoG = CoG, CoY = CoY, response = TRUE)
pred_prob <- predict_lucid(model = early_bin$fit_refit,
                           G = early_bin$inputs$G, Z = early_bin$inputs$Z,
                           CoG = CoG, CoY = CoY, response = FALSE)

cat("response = TRUE  ->", paste(head(pred_lab$pred.y, 8), collapse = " "), "(class labels)\n")
cat("response = FALSE ->", paste(round(head(pred_prob$pred.y, 8), 3), collapse = " "), "(probabilities)\n")

## ----status-table-------------------------------------------------------------
check_obj("early_bin",             "list", "7 early binary")
check_obj("parallel_bin",          "list", "8 parallel binary")
check_obj("serial_ae_bin",         "list", "9 serial all-early binary")
check_obj("serial_mixed_bin",      "list", "10 serial mixed binary")
check_obj("prof_early_bin",        "list", "11 omics profile (early)")
check_obj("prof_par_bin",          "list", "11 omics profile (parallel)")
check_obj("prof_ser_bin",          "list", "11 omics profile (serial)")
check_obj("pred_lab",              "list", "12 predict (labels)")
check_obj("pred_prob",             "list", "12 predict (probabilities)")

status <- do.call(rbind, .reg$rows)
print(status, row.names = FALSE)

cat(sprintf("\n%d of %d registered steps ok; %d not ok\n",
            sum(status$status == "ok"), nrow(status),
            sum(status$status != "ok")))

## -----------------------------------------------------------------------------
sessionInfo()


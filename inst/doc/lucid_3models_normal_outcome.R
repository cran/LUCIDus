## ----setup, message=FALSE, warning=FALSE--------------------------------------
# Keep knitting on error, so that one failing step reports itself and the rest of
# the tutorial still runs. The status table in section 17 records what happened.
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
# which would make the g-computation demonstration in section 12.3 vacuous.
G <- cbind(
  g_causal_1 =  1.0 * (x_true - 0.5) + rnorm(n, sd = 0.8),   # strongest
  g_causal_2 = -0.8 * (x_true - 0.5) + rnorm(n, sd = 0.8),   # moderate, negative
  g_causal_3 =  0.6 * (x_true - 0.5) + rnorm(n, sd = 0.8),   # weakest
  g_noise_1 = rnorm(n), g_noise_2 = rnorm(n), g_noise_3 = rnorm(n),
  g_noise_4 = rnorm(n), g_noise_5 = rnorm(n), g_noise_6 = rnorm(n)
)
G <- as.matrix(scale(G))

causal_exposures <- c("g_causal_1", "g_causal_2", "g_causal_3")

# Exposure penalty used throughout. Section 7.2 shows what this value recovers;
# section 11 sweeps it, and sweeps the omics penalty separately.
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

## ----early-penalized-fit, warning=FALSE---------------------------------------
set.seed(1101)

# Screening model: penalties help identify a parsimonious subset.
early_fit_pen <- estimate_lucid(
  lucid_model = "early",
  G = G,
  Z = Z_early_miss,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  family = "normal",
  K = 2,
  Rho_G = RHO_G,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 15,
  max_tot.itr = 40,
  tol = 1e-2,
  seed = 1101,
  verbose = FALSE
)

# Summary shows model fit, missing profile, selection, and parameter tables.
summary(early_fit_pen)

## ----early-object-structure---------------------------------------------------
names(early_fit_pen)

## ----early-object-fields------------------------------------------------------
cat("log-likelihood:", early_fit_pen$likelihood, "\n")
cat("exposures selected:", sum(get_selected_G(early_fit_pen)),
    "of", length(get_selected_G(early_fit_pen)), "\n")
str(early_fit_pen$em_control, max.level = 1)

## ----early-selected-features--------------------------------------------------
early_selected_G <- names(which(get_selected_G(early_fit_pen)))

cat("exposures retained:", paste(early_selected_G, collapse = ", "), "\n\n")

# Score the selection against the planted truth from section 5: the three
# g_causal_* exposures drive the latent subgroup, the six g_noise_* do not.
data.frame(
  exposure = colnames(G),
  truth    = ifelse(colnames(G) %in% causal_exposures, "causal", "noise"),
  selected = colnames(G) %in% early_selected_G
)

cat("\ncausal exposures kept :", sum(causal_exposures %in% early_selected_G), "of 3\n")
cat("noise  exposures kept :",
    sum(!(early_selected_G %in% causal_exposures)), "of 6\n")

# The omics side is unpenalized in this screening fit, so every feature is
# retained. Section 11.3 penalizes the omics separately and shows why.
cat("omics features retained:", sum(get_selected_Z(early_fit_pen)),
    "of", length(get_selected_Z(early_fit_pen)), "(Rho_Z_Mu = 0 here)\n")

## ----early-selected-refit, warning=FALSE--------------------------------------
set.seed(1102)

# Build selected-only inputs.
early_inputs_refit <- prepare_early_selected_inputs(early_fit_pen, G, Z_early_miss)

# Refit with all penalties set to zero.
early_fit_refit <- refit_selected(
  "early",
  fit_pen = early_fit_pen,
  inputs = early_inputs_refit,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  seed = 1102,
  verbose = TRUE
)

summary(early_fit_refit)

## ----early-bootstrap, warning=FALSE-------------------------------------------
set.seed(1103)

# Bootstrap on zero-penalty refit for CI inference.
early_boot <- boot_lucid(
  G = early_inputs_refit$G,
  Z = early_inputs_refit$Z,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  model = early_fit_refit,
  R = 30,
  conf = 0.90
)

# Print summary with CI columns integrated into parameter tables.
summary(early_fit_refit, boot.se = early_boot)

## ----parallel-penalized-fit, warning=FALSE------------------------------------
set.seed(1201)

parallel_fit_pen <- estimate_lucid(
  lucid_model = "parallel",
  G = G,
  Z = Z_parallel_miss,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  family = "normal",
  K = c(2, 2, 2),
  Rho_G = RHO_G,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 15,
  max_tot.itr = 40,
  tol = 1e-2,
  seed = 1201,
  verbose = FALSE
)

summary(parallel_fit_pen)

## ----parallel-object-structure------------------------------------------------
names(parallel_fit_pen)

## ----parallel-object-fields---------------------------------------------------
cat("log-likelihood:", parallel_fit_pen$likelihood, "\n")
cat("sample size (N):", parallel_fit_pen$N, "\n")
cat("exposures selected per layer:\n")
print(sapply(seq_along(parallel_fit_pen$K), function(i) sum(get_selected_G(parallel_fit_pen, layer = i))))

## ----parallel-selected-features-----------------------------------------------
# Exposure selection union across layers.
parallel_selected_G_union <- names(which(get_selected_G(parallel_fit_pen)))

# Exposure selection per layer.
parallel_selected_G_layer <- lapply(seq_along(parallel_fit_pen$K), function(i) {
  names(which(get_selected_G(parallel_fit_pen, layer = i)))
})

# Omics selection per layer (get_selected_Z already collapses any
# vector/matrix selection encoding to one logical value per feature).
parallel_selected_Z_layer <- lapply(get_selected_Z(parallel_fit_pen), function(s) names(which(s)))

parallel_selected_G_union
parallel_selected_G_layer
parallel_selected_Z_layer

## ----parallel-selected-refit, warning=FALSE-----------------------------------
set.seed(1202)

parallel_inputs_refit <- prepare_parallel_selected_inputs(parallel_fit_pen, G, Z_parallel_miss)
parallel_fit_refit <- refit_selected(
  "parallel",
  fit_pen = parallel_fit_pen,
  inputs = parallel_inputs_refit,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  seed = 1202,
  verbose = TRUE
)

summary(parallel_fit_refit)

## ----parallel-bootstrap, warning=FALSE----------------------------------------
set.seed(1203)

parallel_boot <- boot_lucid(
  G = parallel_inputs_refit$G,
  Z = parallel_inputs_refit$Z,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  model = parallel_fit_refit,
  R = 30,
  conf = 0.90
)

summary(parallel_fit_refit, boot.se = parallel_boot)

## ----serial-all-early-penalized, warning=FALSE--------------------------------
set.seed(1301)

# Serial structure: list of early-stage matrices.
Z_serial_all_early <- list(
  methylome = Z_parallel_miss[[1]],
  transcriptome = Z_parallel_miss[[2]],
  miRNA = Z_parallel_miss[[3]]
)

serial_all_early_pen <- estimate_lucid(
  lucid_model = "serial",
  G = G,
  Z = Z_serial_all_early,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  family = "normal",
  K = list(2, 2, 2),
  Rho_G = RHO_G,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 15,
  max_tot.itr = 40,
  tol = 1e-2,
  seed = 1301,
  verbose = FALSE
)

summary(serial_all_early_pen)
serial_selection_report(serial_all_early_pen)

## ----serial-object-structure--------------------------------------------------
names(serial_all_early_pen)

## ----serial-object-fields-----------------------------------------------------
cat("top-level likelihood (sum over stages):", serial_all_early_pen$likelihood, "\n")
cat("per-stage likelihoods:",
    paste(round(sapply(serial_all_early_pen$submodel, `[[`, "likelihood"), 2),
          collapse = ", "), "\n")
cat("n_stages:", length(serial_all_early_pen$submodel), "\n")

## ----serial-all-early-select, warning=FALSE-----------------------------------
# Top-level select == stage 1's select
identical(serial_all_early_pen$select, serial_all_early_pen$submodel[[1]]$select)

cat("Stage 1 exposures kept:",
    paste(names(which(get_selected_G(serial_all_early_pen))), collapse = ", "), "\n")

# Rho_G is 0 from stage 2 on -- there is nothing there to select
sapply(serial_all_early_pen$submodel, function(sm) sm$Rho$Rho_G)

## ----serial-all-early-refit, warning=FALSE------------------------------------
set.seed(1302)

serial_all_early_inputs_refit <- prepare_serial_selected_inputs(
  serial_all_early_pen, G, Z_serial_all_early
)
serial_all_early_refit <- refit_selected(
  "serial",
  fit_pen = serial_all_early_pen,
  inputs = serial_all_early_inputs_refit,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  seed = 1302,
  verbose = TRUE
)

summary(serial_all_early_refit)

## ----serial-all-early-bootstrap, warning=FALSE--------------------------------
set.seed(1303)

serial_all_early_boot <- boot_lucid(
  G = serial_all_early_inputs_refit$G,
  Z = serial_all_early_inputs_refit$Z,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  model = serial_all_early_refit,
  R = 30,
  conf = 0.90
)

summary(serial_all_early_refit, boot.se = serial_all_early_boot)

## ----serial-mixed-penalized, warning=FALSE------------------------------------
set.seed(1401)

# Nested list signals a parallel submodel at stage 1, followed by early stage 2.
Z_serial_mixed <- list(
  list(
    methylome = Z_parallel_miss[[1]],
    transcriptome = Z_parallel_miss[[2]]
  ),
  miRNA = Z_parallel_miss[[3]]
)

serial_mixed_pen <- estimate_lucid(
  lucid_model = "serial",
  G = G,
  Z = Z_serial_mixed,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  family = "normal",
  K = list(list(2, 2), 2),
  Rho_G = RHO_G,
  Rho_Z_Mu = 0,
  Rho_Z_Cov = 0,
  max_itr = 15,
  max_tot.itr = 40,
  tol = 1e-2,
  seed = 1401,
  verbose = FALSE
)

summary(serial_mixed_pen)
serial_selection_report(serial_mixed_pen)

## ----serial-mixed-refit, warning=FALSE----------------------------------------
set.seed(1402)

serial_mixed_inputs_refit <- prepare_serial_selected_inputs(serial_mixed_pen, G, Z_serial_mixed)
serial_mixed_refit <- refit_selected(
  "serial",
  fit_pen = serial_mixed_pen,
  inputs = serial_mixed_inputs_refit,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  seed = 1402,
  verbose = FALSE
)

summary(serial_mixed_refit)

## ----serial-mixed-bootstrap, warning=FALSE------------------------------------
set.seed(1403)

serial_mixed_boot <- boot_lucid(
  G = serial_mixed_inputs_refit$G,
  Z = serial_mixed_inputs_refit$Z,
  Y = Y,
  CoG = CoG,
  CoY = CoY,
  model = serial_mixed_refit,
  R = 30,
  conf = 0.90
)

summary(serial_mixed_refit, boot.se = serial_mixed_boot)

## ----missing-diagnostics------------------------------------------------------
na_early <- check_na(Z_early_miss, lucid_model = "early")

cat("subjects by missingness code (1 complete, 2 sporadic, 3 listwise):\n")
print(table(na_early$indicator_na))

cat("\nimpute_flag (TRUE only when sporadic subjects exist):",
    na_early$impute_flag, "\n")

## ----missing-pattern----------------------------------------------------------
pat <- analyze_missing_pattern(Z_early_miss)

cat("distinct missingness patterns :", pat$n_patterns, "\n")
cat("complete subjects             :", pat$n_complete, "\n")
cat("overall proportion missing    :", round(pat$total_missing, 4), "\n")

## ----missing-verify-----------------------------------------------------------
fitted_Z <- early_fit_pen$Z

# The listwise subject must still be NA after fitting.
stopifnot(all(is.na(fitted_Z[1, ])))

# The sporadic cells must have been filled.
stopifnot(!anyNA(fitted_Z[2:5, ]))

# Observed values must never have been altered.
obs <- !is.na(Z_early_miss)
stopifnot(isTRUE(all.equal(fitted_Z[obs], Z_early_miss[obs])))

cat("verified: listwise row still NA, sporadic cells filled,",
    "observed values unchanged\n")

## ----missing-serial-note------------------------------------------------------
cat("early  fit$Z still has NAs? ", anyNA(early_fit_pen$Z), "\n")
cat("serial fit$Z still has NAs? ", anyNA(serial_all_early_pen$Z[[1]]),
    "  (imputed data lives on submodel[[i]]$Z)\n")

## ----safe-impute--------------------------------------------------------------
Z_mean_filled <- safe_impute(Z_early_miss, method = "mean")

qual <- check_imputation_quality(Z_early_miss, Z_mean_filled)
cat("mean-filled  : valid =", qual$is_valid,
    " mean shift =", round(qual$mean_diff, 3),
    " sd ratio =", round(qual$sd_ratio, 3), "\n")

## ----penalty-sweep-z, warning=FALSE-------------------------------------------
sweep_rho_z <- function(rz) {
  invisible(capture.output(
    f <- estimate_lucid(
      lucid_model = "early", G = G, Z = Z_early_miss, Y = Y, CoG = CoG, CoY = CoY,
      family = "normal", K = 2, init_omic.data.model = "EEV",
      Rho_G = 0, Rho_Z_Mu = rz, Rho_Z_Cov = if (rz > 0) 0.02 else 0,
      max_itr = 30, max_tot.itr = 90, tol = 1e-2, seed = 1101
    )
  ))
  kept <- which(get_selected_Z(f))
  acc <- mean(get_cluster_assignment(f) == x_true + 1); acc <- max(acc, 1 - acc)
  data.frame(
    Rho_Z_Mu       = rz,
    features_kept  = length(kept),
    signal_kept    = sum(kept %in% signal_cols_early),
    converged      = isTRUE(f$em_control$converged),
    subgroup_recovery = round(acc, 2)
  )
}

do.call(rbind, lapply(c(0, 2, 10), sweep_rho_z))

## ----predict-training---------------------------------------------------------
# CoG and CoY must be supplied here exactly as they were at fitting time: the
# G -> X design includes the CoG columns, so omitting them makes the design
# matrix too narrow for the fitted coefficients.
# lucid_model is auto-detected from early_fit_refit's own class.
pred_train <- predict_lucid(
  model = early_fit_refit,
  G = early_inputs_refit$G, Z = early_inputs_refit$Z, Y = Y,
  CoG = CoG, CoY = CoY
)

# Predicting on the training data reproduces the fit's own assignment.
stopifnot(identical(
  as.numeric(pred_train$pred.x),
  as.numeric(get_cluster_assignment(early_fit_refit))
))

cat("cluster sizes:\n"); print(table(pred_train$pred.x))

## ----predict-legality---------------------------------------------------------
legal <- function(expr) {
  r <- try(suppressWarnings(suppressMessages(invisible(capture.output(force(expr))))),
           silent = TRUE)
  if (inherits(r, "try-error")) "error" else "works"
}

Gp <- early_inputs_refit$G
Zp <- early_inputs_refit$Z
call_p <- function(...) predict_lucid(model = early_fit_refit,
                                      G = Gp, CoG = CoG, CoY = CoY, ...)

data.frame(
  Z              = c("supplied", "supplied", "omitted", "omitted", "omitted"),
  Y              = c("supplied", "omitted",  "supplied", "omitted", "omitted"),
  g_computation  = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  result         = c(
    legal(call_p(Z = Zp, Y = Y)),
    legal(call_p(Z = Zp)),
    legal(call_p(Z = NULL, Y = Y)),
    legal(call_p(Z = NULL)),
    legal(call_p(Z = NULL, g_computation = TRUE))
  )
)

## ----predict-missing-z-message------------------------------------------------
cat(tryCatch(
  call_p(Z = NULL, Y = Y),
  error = function(e) conditionMessage(e)
))

## ----predict-all-models, warning=FALSE----------------------------------------
# pred.x nests: a vector for early, a list by layer for parallel, and for serial
# a list by stage whose elements are themselves lists when that stage is a
# parallel submodel. Flatten to the leaves so each cluster variable is counted on
# its own -- pooling a parallel stage's layers would report twice as many
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

outcome_summary <- function(pred_y, label) {
  v <- as.numeric(unlist(pred_y))
  data.frame(model = label, n = length(v),
             mean = round(mean(v), 3), sd = round(stats::sd(v), 3),
             min = round(min(v), 3), median = round(stats::median(v), 3),
             max = round(max(v), 3), row.names = NULL)
}

pred_early_m <- predict_lucid(
  model = early_fit_refit,
  G = early_inputs_refit$G, Z = early_inputs_refit$Z, Y = Y,
  CoG = CoG, CoY = CoY
)
pred_parallel_m <- predict_lucid(
  model = parallel_fit_refit,
  G = parallel_inputs_refit$G, Z = parallel_inputs_refit$Z, Y = Y,
  CoG = CoG, CoY = CoY
)
pred_serial_m <- predict_lucid(
  model = serial_all_early_refit,
  G = serial_all_early_inputs_refit$G, Z = serial_all_early_inputs_refit$Z, Y = Y,
  CoG = CoG, CoY = CoY
)

rbind(
  cluster_sizes(pred_early_m$pred.x,    "early"),
  cluster_sizes(pred_parallel_m$pred.x, "parallel (per layer)"),
  cluster_sizes(pred_serial_m$pred.x,   "serial (per stage)")
)

## ----predict-all-models-y-----------------------------------------------------
rbind(
  outcome_summary(pred_early_m$pred.y,    "early"),
  outcome_summary(pred_parallel_m$pred.y, "parallel"),
  outcome_summary(pred_serial_m$pred.y,   "serial"),
  outcome_summary(Y,                      "observed outcome")
)

## ----predict-supervised-------------------------------------------------------
pred_unsup <- predict_lucid(
  model = early_fit_refit,
  G = early_inputs_refit$G, Z = early_inputs_refit$Z,
  CoG = CoG, CoY = CoY
)

cat("agreement between supervised and unsupervised assignment:",
    round(mean(pred_train$pred.x == pred_unsup$pred.x), 3), "\n")

## ----g-computation------------------------------------------------------------
G_ref <- early_inputs_refit$G

# Shift every retained causal exposure by 1.5 SD, each in the direction of its
# own effect, i.e. "what if this whole exposure profile were less favourable?"
kept_causal <- intersect(colnames(G_ref), causal_exposures)
direction <- c(g_causal_1 = 1, g_causal_2 = -1, g_causal_3 = 1)

shift_profile <- function(Gx, delta) {
  for (nm in kept_causal) Gx[, nm] <- Gx[, nm] + delta * direction[[nm]]
  Gx
}

gc_hi <- predict_lucid(model = early_fit_refit,
                       G = shift_profile(G_ref,  1.5), Z = NULL,
                       CoG = CoG, CoY = CoY, g_computation = TRUE)
gc_lo <- predict_lucid(model = early_fit_refit,
                       G = shift_profile(G_ref, -1.5), Z = NULL,
                       CoG = CoG, CoY = CoY, g_computation = TRUE)

data.frame(
  scenario          = c("profile +1.5 SD", "profile -1.5 SD"),
  mean_pred_outcome = c(mean(gc_hi$pred.y), mean(gc_lo$pred.y)),
  prop_in_cluster_2 = c(mean(gc_hi$inclusion.p[, 2]), mean(gc_lo$inclusion.p[, 2]))
)

## ----g-computation-decompose--------------------------------------------------
# early_gamma_levels() is the internal accessor predict_lucid() itself uses to
# turn the fitted outcome parameters into one level per cluster -- not part of
# the public API, called here (via :::) only to show the arithmetic behind the
# g-computation contrast above. Reading res_Gamma$beta directly would give the
# stored parameterization, which is not the same thing.
levels_by_cluster <- LUCIDus:::early_gamma_levels(early_fit_refit$res_Gamma, early_fit_refit$K)
shift_in_membership <- mean(gc_hi$inclusion.p[, 2]) - mean(gc_lo$inclusion.p[, 2])

cat("cluster outcome levels     :", round(levels_by_cluster, 3), "\n")
cat("gap between clusters       :", round(diff(levels_by_cluster), 3), "\n")
cat("shift in cluster-2 share   :", round(shift_in_membership, 3), "\n")
cat("implied outcome contrast   :",
    round(shift_in_membership * diff(levels_by_cluster), 4), "\n")
cat("observed outcome contrast  :",
    round(mean(gc_hi$pred.y) - mean(gc_lo$pred.y), 4), "\n")

## ----g-computation-pred-z-----------------------------------------------------
cat("predicted omics means differ between scenarios by (first 5 features):\n")
print(round(head(colMeans(gc_hi$pred.z) - colMeans(gc_lo$pred.z), 5), 4))

## ----lucid-wrapper, warning=FALSE---------------------------------------------
set.seed(1501)
fit_tuned <- lucid(
  G = G, Z = Z_early_miss, Y = Y, CoG = CoG, CoY = CoY,
  lucid_model = "early", family = "normal",
  K = 2:3,
  Rho_G = c(0, 0.02),
  init_omic.data.model = NULL,
  max_itr = 15, max_tot.itr = 40, tol = 1e-2
)

cat("selected K :", fit_tuned$K, "\n")
cat("BIC        :", round(summary(fit_tuned, auto_print = FALSE)$BIC, 2), "\n")

## ----lucid-selection----------------------------------------------------------
if (!is.null(fit_tuned$selection)) {
  cat("exposures retained:",
      paste(fit_tuned$selection$Gnames[fit_tuned$selection$selectG], collapse = ", "), "\n")
  cat("exposures dropped :",
      paste(fit_tuned$selection$Gnames[!fit_tuned$selection$selectG], collapse = ", "), "\n")
} else {
  cat("no variable was deselected at this penalty grid\n")
}

## ----sankey-------------------------------------------------------------------
plot(early_fit_refit)

## ----profile-early-heat, fig.width = 6.5, fig.height = 4.5, warning=FALSE-----
prof_early <- plot_cluster_omic_profile(early_fit_refit, top_n = 10)
prof_early[[1]]

## ----profile-importance-compare, warning=FALSE--------------------------------
rank_by <- function(measure) {
  d <- attr(plot_cluster_omic_profile(early_fit_refit, top_n = 5,
                                      importance = measure)[[1]], "profile_data")
  as.character(unique(d$feature)[order(-unique(d[, c("feature","score")])$score)])
}

data.frame(
  rank        = 1:5,
  separation  = rank_by("separation"),
  range       = rank_by("range")
)

## ----profile-early-bar, fig.width = 6.5, fig.height = 4.5, warning=FALSE------
plot_cluster_omic_profile(early_fit_refit, type = "bar", top_n = 10)[[1]]

## ----profile-parallel, fig.width = 6.5, fig.height = 4.5, warning=FALSE-------
prof_par <- plot_cluster_omic_profile(
  parallel_fit_refit,
  layer_names = c("methylome", "transcriptome", "miRNA"),
  top_n = 8
)
names(prof_par)
for (nm in names(prof_par)) print(prof_par[[nm]])

## ----profile-parallel-bar, fig.width = 6.5, fig.height = 4.5, warning=FALSE----
plot_cluster_omic_profile(parallel_fit_refit, type = "bar",
                          layer_names = c("methylome", "transcriptome", "miRNA"),
                          top_n = 8)[["transcriptome"]]

## ----profile-serial, fig.width = 6.5, fig.height = 4.5, warning=FALSE---------
prof_ser <- plot_cluster_omic_profile(serial_all_early_refit, top_n = 8)
names(prof_ser)
for (nm in names(prof_ser)) print(prof_ser[[nm]])

## ----profile-serial-bar, fig.width = 6.5, fig.height = 4.5, warning=FALSE-----
plot_cluster_omic_profile(serial_all_early_refit, type = "bar", top_n = 8)[[2]]

## ----profile-data, warning=FALSE----------------------------------------------
pd <- attr(prof_early[[1]], "profile_data")
head(pd[order(-pd$score), c("feature", "cluster", "mean", "sd", "score")], 6)

## ----status-table-------------------------------------------------------------
check_obj("early_fit_refit",       "early_lucid",    "7 early normal")
check_obj("parallel_fit_refit",    "lucid_parallel", "8 parallel normal")
check_obj("serial_all_early_refit","lucid_serial",   "9 serial all-early normal")
check_obj("serial_mixed_refit",    "lucid_serial",   "10 serial mixed normal")
check_obj("na_early",              "list",           "11 missing diagnostics")
check_obj("pred_train",            "list",           "12 prediction")
check_obj("gc_hi",                 "list",           "12.3 g-computation")
check_obj("fit_tuned",             "early_lucid",    "13 lucid() wrapper")
check_obj("prof_early",            "list",           "15 omics profile (early)")
check_obj("prof_par",              "list",           "15.3 omics profile (parallel)")
check_obj("prof_ser",              "list",           "15.3 omics profile (serial)")

status <- do.call(rbind, .reg$rows)
print(status, row.names = FALSE)

cat(sprintf("\n%d of %d registered steps ok; %d not ok\n",
            sum(status$status == "ok"), nrow(status),
            sum(status$status != "ok")))

## -----------------------------------------------------------------------------
sessionInfo()


## ----setup, message=FALSE, warning=FALSE--------------------------------------
required_pkgs <- c("plotly")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install required packages first: ", paste(missing_pkgs, collapse = ", "))
}

library(LUCIDus)
library(plotly)  # for %>% and plot_ly(), used directly by this vignette's own Sankey code

data(sim_data)

## ----simulate-omics-----------------------------------------------------------
set.seed(1008)

n <- 200

G <- sim_data$G[seq_len(n), , drop = FALSE]
colnames(G) <- paste0("exposure", seq_len(ncol(G)))

Y <- as.numeric(sim_data$Y_normal[seq_len(n), 1])
CoY <- sim_data$Covariate[seq_len(n), , drop = FALSE]
colnames(CoY) <- paste0("covariate", seq_len(ncol(CoY)))

latent_signal <- ifelse(sim_data$X[seq_len(n)] == 1, -1, 1)
g_anchor <- as.numeric(scale(G[, 1]))

simulate_omics_block <- function(feature_names, signal_n, latent_signal, g_anchor,
                                 cluster_scale, noise_sd, offset = 0) {
  n_obs <- length(latent_signal)
  p <- length(feature_names)
  block <- matrix(NA_real_, nrow = n_obs, ncol = p)

  signal_loadings <- seq(cluster_scale, 0.8, length.out = signal_n)

  for (j in seq_len(p)) {
    direction <- if (j %% 2L == 0L) 1 else -1
    cluster_loading <- if (j <= signal_n) {
      direction * signal_loadings[j]
    } else {
      runif(1, -0.20, 0.20)
    }
    exposure_loading <- if (j <= signal_n) {
      direction * runif(1, 0.15, 0.35)
    } else {
      runif(1, -0.10, 0.10)
    }

    block[, j] <- offset +
      cluster_loading * latent_signal +
      exposure_loading * g_anchor +
      rnorm(n_obs, sd = noise_sd)
  }

  block <- scale(block)
  colnames(block) <- feature_names
  block
}

n_microbiota <- 18
n_hilic <- 12
n_c18 <- 8

microbiota_names <- paste0("X", sprintf("%07d", seq_len(n_microbiota)))
hilic_names <- paste0("hilic", sprintf("%04d", seq_len(n_hilic)))
c18_names <- paste0("c18", sprintf("%04d", seq_len(n_c18)))

microbiota <- simulate_omics_block(
  feature_names = microbiota_names,
  signal_n = 6,
  latent_signal = latent_signal,
  g_anchor = g_anchor,
  cluster_scale = 1.80,
  noise_sd = 0.85
)

hilic <- simulate_omics_block(
  feature_names = hilic_names,
  signal_n = 5,
  latent_signal = latent_signal,
  g_anchor = g_anchor,
  cluster_scale = 1.55,
  noise_sd = 0.90,
  offset = 0.15
)

c18 <- simulate_omics_block(
  feature_names = c18_names,
  signal_n = 4,
  latent_signal = latent_signal,
  g_anchor = g_anchor,
  cluster_scale = 1.35,
  noise_sd = 0.95,
  offset = -0.10
)

metabolites <- cbind(hilic, c18)
Z <- cbind(microbiota, metabolites)

omics_dimensions <- data.frame(
  block = c("microbiota", "hilic", "c18", "metabolites_total", "combined_Z"),
  n_features = c(ncol(microbiota), ncol(hilic), ncol(c18), ncol(metabolites), ncol(Z))
)
omics_dimensions

## ----fit-early, message=FALSE-------------------------------------------------
fit_early <- estimate_lucid(
  G = G,
  Z = Z,
  Y = Y,
  lucid_model = "early",
  CoY = CoY,
  family = "normal",
  K = 2,
  seed = 1008,
  max_itr = 20,
  max_tot.itr = 60,
  tol = 1e-2
)

class(fit_early)
dim(fit_early$res_Mu)

## ----sankey-function----------------------------------------------------------
detect_omics_family <- function(x) {
  out <- rep("exposure", length(x))
  out[grepl("^Latent Cluster", x)] <- "lc"
  out[grepl("^X\\d{7}$", x)] <- "microbiota"
  out[grepl("^(hilic|c18)", x, ignore.case = TRUE)] <- "metabolite"
  out
}

select_top_omics_balanced <- function(mu_mat, z_names, cluster_pair, top_n_each = NULL) {
  if (is.null(top_n_each)) {
    return(z_names)
  }

  if (!is.numeric(top_n_each) || length(top_n_each) != 1 || is.na(top_n_each) || top_n_each < 1) {
    stop("`top_n_each` must be a single positive integer or NULL.")
  }

  top_n_each <- as.integer(top_n_each)
  c1 <- as.integer(cluster_pair[1])
  c2 <- as.integer(cluster_pair[2])
  score <- abs(mu_mat[c1, ] - mu_mat[c2, ])
  family <- detect_omics_family(z_names)

  ordered_keep <- character(0)
  for (fam in c("microbiota", "metabolite")) {
    fam_idx <- which(family == fam)
    if (!length(fam_idx)) {
      next
    }
    fam_idx <- fam_idx[order(score[fam_idx], decreasing = TRUE)]
    ordered_keep <- c(
      ordered_keep,
      z_names[fam_idx[seq_len(min(top_n_each, length(fam_idx)))]]
    )
  }

  ordered_keep
}

make_link_color <- function(target_group, is_positive) {
  pos_cols <- c(
    lc = "#8E9BAE",
    microbiota = "#1B7F5A",
    metabolite = "#C46A1B",
    outcome = "#B23A48",
    exposure = "#5E6F86"
  )
  neg_cols <- c(
    lc = "#DCE3EA",
    microbiota = "#BFE2D6",
    metabolite = "#F2D4B5",
    outcome = "#E9C0C6",
    exposure = "#CCD6E2"
  )

  ifelse(is_positive, pos_cols[target_group], neg_cols[target_group])
}

sankey_early_integration <- function(
    lucid_fit1,
    top_n_each = 5,
    cluster_pair = c(1, 2),
    text_size = 13,
    drop_covariates = TRUE
) {
  if (!inherits(lucid_fit1, "early_lucid")) {
    stop("`lucid_fit1` must inherit from class 'early_lucid'.")
  }

  needed <- c("K", "var.names", "res_Beta", "res_Mu", "res_Gamma")
  miss <- needed[!needed %in% names(lucid_fit1)]
  if (length(miss) > 0) {
    stop("Missing fields in `lucid_fit1`: ", paste(miss, collapse = ", "))
  }

  K <- lucid_fit1$K
  var.names <- lucid_fit1$var.names
  Gnames <- var.names$Gnames
  Znames <- var.names$Znames
  Yname <- var.names$Ynames

  beta_mat <- as.matrix(lucid_fit1$res_Beta)
  mu_mat <- as.matrix(lucid_fit1$res_Mu)

  if (nrow(beta_mat) != K) stop("Unexpected `res_Beta` shape.")
  if (nrow(mu_mat) != K) stop("Unexpected `res_Mu` shape.")

  if (length(cluster_pair) != 2) stop("`cluster_pair` must be length-2.")
  c1 <- as.integer(cluster_pair[1])
  c2 <- as.integer(cluster_pair[2])
  if (c1 < 1 || c1 > K || c2 < 1 || c2 > K || c1 == c2) {
    stop("`cluster_pair` must reference two distinct cluster indices in 1..K.")
  }

  keep_z <- select_top_omics_balanced(
    mu_mat = mu_mat,
    z_names = Znames,
    cluster_pair = cluster_pair,
    top_n_each = top_n_each
  )

  valueGtoX <- as.vector(t(beta_mat[, -1, drop = FALSE]))

  mu_use <- mu_mat[, Znames %in% keep_z, drop = FALSE]
  mu_use <- mu_use[, keep_z, drop = FALSE]
  z_use <- colnames(mu_use)
  valueXtoZ <- as.vector(t(mu_use))

  gamma_beta <- if (!is.null(lucid_fit1$res_Gamma$beta)) {
    as.numeric(lucid_fit1$res_Gamma$beta)
  } else {
    as.numeric(lucid_fit1$res_Gamma$fit$coefficients)
  }
  if (length(gamma_beta) < K) stop("Could not extract cluster->outcome effects.")

  valueXtoY <- gamma_beta[seq_len(K)]
  valueXtoY[1] <- 0

  GtoX <- data.frame(
    source = rep(Gnames, K),
    target = paste0("Latent Cluster", rep(seq_len(K), each = length(Gnames))),
    value = abs(valueGtoX),
    group = valueGtoX > 0,
    stringsAsFactors = FALSE
  )

  XtoZ <- data.frame(
    source = paste0("Latent Cluster", rep(seq_len(K), each = length(z_use))),
    target = rep(z_use, K),
    value = abs(valueXtoZ),
    group = valueXtoZ > 0,
    stringsAsFactors = FALSE
  )

  XtoY <- data.frame(
    source = paste0("Latent Cluster", seq_len(K)),
    target = rep(Yname, K),
    value = abs(valueXtoY),
    group = valueXtoY > 0,
    stringsAsFactors = FALSE
  )

  links <- rbind(GtoX, XtoZ, XtoY)

  if (drop_covariates) {
    links <- links[!grepl("^covariate", links$source, ignore.case = TRUE), , drop = FALSE]
  }

  nodes <- data.frame(name = unique(c(links$source, links$target)), stringsAsFactors = FALSE)
  nodes$group <- detect_omics_family(nodes$name)
  nodes$group[nodes$name %in% Yname] <- "outcome"

  group_order <- c("exposure", "lc", "microbiota", "metabolite", "outcome")
  nodes$group <- factor(nodes$group, levels = group_order)
  nodes <- nodes[order(nodes$group, nodes$name), , drop = FALSE]
  nodes$group <- as.character(nodes$group)

  node_color_map <- c(
    exposure = "#4C6A92",
    lc = "#B9C6D3",
    microbiota = "#2A9D6F",
    metabolite = "#D98C3F",
    outcome = "#C44536"
  )
  nodes$color <- unname(node_color_map[nodes$group])

  links$IDsource <- match(links$source, nodes$name) - 1L
  links$IDtarget <- match(links$target, nodes$name) - 1L

  target_group <- nodes$group[match(links$target, nodes$name)]
  target_group[grepl("^Latent Cluster", links$target)] <- "lc"
  links$link_color <- make_link_color(target_group, links$group)

  x_map <- c(
    exposure = 0.02,
    lc = 0.34,
    microbiota = 0.68,
    metabolite = 0.80,
    outcome = 0.97
  )
  nodes$x <- unname(x_map[nodes$group])

  nodes$y <- NA_real_
  for (g in group_order) {
    idx <- which(nodes$group == g)
    if (!length(idx)) {
      next
    }
    nodes$y[idx] <- if (length(idx) == 1) 0.5 else seq(0.03, 0.97, length.out = length(idx))
  }

  plotly::plot_ly(
    type = "sankey",
    orientation = "h",
    node = list(
      label = nodes$name,
      color = nodes$color,
      pad = 16,
      thickness = 18,
      line = list(color = "black", width = 0.35),
      x = nodes$x,
      y = nodes$y
    ),
    link = list(
      source = links$IDsource,
      target = links$IDtarget,
      value = pmax(links$value, 1e-12),
      color = links$link_color
    )
  ) %>%
    plotly::layout(
      font = list(size = text_size),
      margin = list(l = 20, r = 20, t = 30, b = 20)
    )
}

## ----selected-features--------------------------------------------------------
selected_top <- select_top_omics_balanced(
  mu_mat = fit_early$res_Mu,
  z_names = fit_early$var.names$Znames,
  cluster_pair = c(1, 2),
  top_n_each = 5
)

data.frame(
  feature = selected_top,
  family = detect_omics_family(selected_top),
  stringsAsFactors = FALSE
)

## ----top-omics-features-------------------------------------------------------
top_general <- get_top_omics_features(fit_early, top_n = 10)
top_general

## ----plot-all, fig.width=13, fig.height=9-------------------------------------
fig_all <- sankey_early_integration(
  lucid_fit1 = fit_early,
  top_n_each = NULL,
  cluster_pair = c(1, 2),
  text_size = 11
)
fig_all

## ----plot-top, fig.width=13, fig.height=8-------------------------------------
fig_top_balanced <- sankey_early_integration(
  lucid_fit1 = fit_early,
  top_n_each = 5,
  cluster_pair = c(1, 2),
  text_size = 13
)
fig_top_balanced


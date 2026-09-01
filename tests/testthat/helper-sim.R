# ---------------------------------------------------------------------------
# Canonical data generator for LUCID validation.
#
# Samples directly from the DAG of Figure 1 (RJ-2024-012) with a known Theta,
# and returns that Theta alongside the data so every test can compare against
# ground truth rather than against the package's own output.
#
# Missingness is imposed AFTER generation and the complete omics matrix is
# retained as `Z_complete`, which is what makes imputation accuracy measurable.
# ---------------------------------------------------------------------------

# softmax rows of a matrix
sim_softmax <- function(eta) {
  t(apply(eta, 1, function(v) {
    v <- v - max(v)
    e <- exp(v)
    e / sum(e)
  }))
}

# draw one categorical value per row of a probability matrix
sim_rcat <- function(prob) {
  vapply(seq_len(nrow(prob)), function(i) {
    sample.int(ncol(prob), size = 1L, prob = prob[i, ])
  }, integer(1))
}

# multivariate normal draws (chol-based, no external dependency)
sim_rmvnorm <- function(n, mean, sigma) {
  p <- length(mean)
  ch <- chol(sigma)
  matrix(stats::rnorm(n * p), nrow = n) %*% ch +
    matrix(rep(mean, each = n), nrow = n)
}

# Default well-separated truth for a K-cluster, M-feature omics block.
sim_default_mu <- function(K, M, sep = 1.5) {
  mu <- matrix(0, nrow = K, ncol = M)
  for (k in seq_len(K)) mu[k, ] <- sep * ((k - 1) - (K - 1) / 2) * rep(c(1, -1), length.out = M)
  mu
}

sim_default_sigma <- function(K, M, sd = 1) {
  sig <- array(0, dim = c(M, M, K))
  for (k in seq_len(K)) sig[, , k] <- diag(sd^2, M)
  sig
}

# Default beta: K x (1 + P), row 1 is the zero reference row.
sim_default_beta <- function(K, P, effect = 1) {
  b <- matrix(0, nrow = K, ncol = P + 1L)
  for (k in seq.int(2, K)) {
    b[k, 1] <- 0
    b[k, 2] <- effect * (k - 1)          # first exposure carries signal
    if (P > 1) b[k, seq.int(3, P + 1L)] <- 0   # remaining exposures are null
  }
  b
}

#' Simulate data from the LUCID DAG.
#'
#' @param model "early", "parallel" or "serial"
#' @param N sample size
#' @param K number of clusters; a vector (one per layer) for parallel/serial
#' @param P number of exposures
#' @param M number of omics features per layer
#' @param n_layer number of omics layers (parallel/serial)
#' @param family "normal" or "binary"
#' @param missing "none", "listwise", "sporadic" or "mixed"
#' @param miss_ratio proportion of rows (listwise) or cells (sporadic) missing
#' @param mechanism "MCAR" or "MAR"; MAR makes missingness depend on G and Y
#' @return list with G, Z, Y, Z_complete and `truth`
sim_lucid <- function(model = c("early", "parallel", "serial"),
                      N = 400, K = 2, P = 2, M = 4, n_layer = 2,
                      family = c("normal", "binary"),
                      beta = NULL, mu = NULL, sigma = NULL,
                      gamma = NULL, gamma_sd = NULL,
                      n_CoG = 0, n_CoY = 0,
                      coef_CoG = 0.5, coef_CoY = 0.4,
                      sep = 1.5, effect = 1,
                      missing = c("none", "listwise", "sporadic", "mixed"),
                      miss_ratio = 0.2,
                      mechanism = c("MCAR", "MAR"),
                      seed = 1L) {
  model <- match.arg(model)
  family <- match.arg(family)
  missing <- match.arg(missing)
  mechanism <- match.arg(mechanism)
  set.seed(seed)

  if (model == "early") n_layer <- 1L
  K <- rep_len(K, n_layer)

  G <- matrix(stats::rnorm(N * P), nrow = N)
  colnames(G) <- paste0("G", seq_len(P))
  CoG <- NULL
  if (n_CoG > 0) {
    CoG <- matrix(stats::rnorm(N * n_CoG), nrow = N)
    colnames(CoG) <- paste0("CoG", seq_len(n_CoG))
  }
  CoY <- NULL
  if (n_CoY > 0) {
    CoY <- matrix(stats::rnorm(N * n_CoY), nrow = N)
    colnames(CoY) <- paste0("CoY", seq_len(n_CoY))
  }

  P_full <- P + n_CoG
  G_full <- if (is.null(CoG)) G else cbind(G, CoG)

  if (is.null(beta)) {
    beta <- lapply(seq_len(n_layer), function(a) {
      b <- sim_default_beta(K[a], P_full, effect = effect)
      if (n_CoG > 0) {
        for (k in seq.int(2, K[a])) b[k, P + 1L + seq_len(n_CoG)] <- coef_CoG
      }
      b
    })
  } else if (!is.list(beta)) beta <- list(beta)

  if (is.null(mu))    mu    <- lapply(seq_len(n_layer), function(a) sim_default_mu(K[a], M, sep))
  else if (!is.list(mu)) mu <- list(mu)
  if (is.null(sigma)) sigma <- lapply(seq_len(n_layer), function(a) sim_default_sigma(K[a], M))
  else if (!is.list(sigma)) sigma <- list(sigma)

  # ---- latent clusters and omics -----------------------------------------
  X <- vector("list", n_layer)
  Z <- vector("list", n_layer)
  delta <- vector("list", max(0L, n_layer - 1L))
  prev <- NULL
  for (a in seq_len(n_layer)) {
    if (model == "serial" && a > 1L) {
      # upstream cluster indicator plays the role of "G" for this stage
      Xprev <- prev
      design <- cbind(1, stats::model.matrix(~ factor(Xprev, levels = seq_len(K[a - 1L])))[, -1, drop = FALSE])
      b <- matrix(0, nrow = K[a], ncol = ncol(design))
      for (k in seq.int(2, K[a])) b[k, -1] <- effect * (k - 1)
      # Record the coefficients actually used for this transition. `beta[[a]]`
      # is the G-based matrix, which a serial stage a > 1 never uses -- it does
      # not even have the right number of columns -- so recovery tests need this
      # separately or they compare against a parameter that was never applied.
      delta[[a - 1L]] <- b
      eta <- design %*% t(b)
    } else {
      eta <- cbind(1, G_full) %*% t(beta[[a]])
    }
    prob <- sim_softmax(eta)
    X[[a]] <- sim_rcat(prob)
    prev <- X[[a]]

    Za <- matrix(0, nrow = N, ncol = M)
    for (k in seq_len(K[a])) {
      idx <- which(X[[a]] == k)
      if (length(idx)) Za[idx, ] <- sim_rmvnorm(length(idx), mu[[a]][k, ], sigma[[a]][, , k])
    }
    colnames(Za) <- paste0("Z", a, "_", seq_len(M))
    Z[[a]] <- Za
  }

  # ---- outcome ------------------------------------------------------------
  if (is.null(gamma)) {
    gamma <- lapply(seq_len(n_layer), function(a) seq_len(K[a]) - 1)   # 0, 1, 2, ...
  } else if (!is.list(gamma)) gamma <- list(gamma)
  if (is.null(gamma_sd)) gamma_sd <- lapply(seq_len(n_layer), function(a) rep(1, K[a]))
  else if (!is.list(gamma_sd)) gamma_sd <- list(gamma_sd)

  if (model == "serial") {
    eta_y <- gamma[[n_layer]][X[[n_layer]]]
  } else {
    eta_y <- rowSums(vapply(seq_len(n_layer), function(a) gamma[[a]][X[[a]]], numeric(N)))
  }
  if (n_CoY > 0) eta_y <- eta_y + as.numeric(CoY %*% rep(coef_CoY, n_CoY))

  if (family == "normal") {
    sd_i <- if (model == "early") gamma_sd[[1]][X[[1]]] else rep(1, N)
    Y <- stats::rnorm(N, eta_y, sd_i)
  } else {
    Y <- stats::rbinom(N, 1, stats::plogis(eta_y))
  }
  Y <- as.matrix(Y); colnames(Y) <- "outcome"

  Z_complete <- lapply(Z, function(x) x)

  # ---- missingness --------------------------------------------------------
  miss_index <- vector("list", n_layer)
  if (missing != "none") {
    # row-selection probability: MCAR is uniform, MAR depends on G and Y
    if (mechanism == "MAR") {
      lp <- 0.9 * G[, 1] + 0.9 * as.numeric(scale(as.numeric(Y)))
      w <- stats::plogis(lp)
    } else {
      w <- rep(1, N)
    }
    for (a in seq_len(n_layer)) {
      Za <- Z[[a]]
      if (missing %in% c("listwise", "mixed")) {
        n_lw <- max(1L, round(miss_ratio * N))
        rows <- sample.int(N, n_lw, prob = w / sum(w))
        Za[rows, ] <- NA
      }
      if (missing %in% c("sporadic", "mixed")) {
        pool <- which(!is.na(Za[, 1]))          # do not touch listwise rows
        n_sp <- max(1L, round(miss_ratio * length(pool) * M * 0.5))
        n_sp <- min(n_sp, length(pool) * (M - 1L))
        cells <- cbind(sample(pool, n_sp, replace = TRUE),
                       sample.int(M, n_sp, replace = TRUE))
        for (i in seq_len(nrow(cells))) {
          # never blank out an entire row via the sporadic path
          if (sum(is.na(Za[cells[i, 1], ])) < M - 1L) Za[cells[i, 1], cells[i, 2]] <- NA
        }
      }
      Z[[a]] <- Za
      miss_index[[a]] <- is.na(Za)
    }
  } else {
    for (a in seq_len(n_layer)) miss_index[[a]] <- matrix(FALSE, N, M)
  }

  out_Z          <- if (model == "early") Z[[1]] else Z
  out_Z_complete <- if (model == "early") Z_complete[[1]] else Z_complete
  out_miss       <- if (model == "early") miss_index[[1]] else miss_index

  list(
    G = G, Z = out_Z, Y = Y, CoG = CoG, CoY = CoY,
    Z_complete = out_Z_complete,
    miss_index = out_miss,
    truth = list(
      model = model, family = family, K = K, N = N, P = P, M = M,
      n_layer = n_layer, beta = beta, mu = mu, sigma = sigma,
      # `delta` is populated for serial only: element i holds the coefficients
      # generating stage i+1 from stage i. Empty list for early/parallel.
      delta = if (model == "serial") delta else list(),
      gamma = gamma, gamma_sd = gamma_sd, X = X,
      n_CoG = n_CoG, n_CoY = n_CoY,
      coef_CoG = coef_CoG, coef_CoY = coef_CoY,
      mechanism = mechanism, missing = missing, miss_ratio = miss_ratio
    )
  )
}

# ---------------------------------------------------------------------------
# Serial data with an explicit per-stage topology.
#
# `sim_lucid("serial", ...)` takes a flat K vector, so every stage it builds is
# single-layer -- an all-early chain. It cannot express a serial model with a
# *parallel* stage inside, which is a topology the package fully supports:
# R/EM_all.R:461 decides a stage is parallel exactly when `K[[i]]` is a list.
#
# `topology` mirrors the `K` argument the package expects:
#   list(2, list(2, 2), 3)  ->  early(K=2) -> parallel(K=2,2) -> early(K=3)
# and `Z` is emitted nested to match, as normalize_serial_block() requires.
#
# Ground truth retained per stage: cluster assignment, mu, sigma, and the
# transition coefficients actually used.
# ---------------------------------------------------------------------------
sim_lucid_serial_topology <- function(topology,
                                      N = 400, P = 2, M = 4,
                                      family = c("normal", "binary"),
                                      n_CoG = 0, n_CoY = 0,
                                      coef_CoG = 0.5, coef_CoY = 0.4,
                                      sep = 1.5, effect = 1,
                                      missing = c("none", "listwise", "sporadic", "mixed"),
                                      miss_ratio = 0.2,
                                      mechanism = c("MCAR", "MAR"),
                                      seed = 1L) {
  family <- match.arg(family)
  missing <- match.arg(missing)
  mechanism <- match.arg(mechanism)
  set.seed(seed)

  n_stage <- length(topology)
  # per stage: the vector of layer cluster counts, and whether it is parallel
  stage_K <- lapply(topology, function(k) as.numeric(unlist(k, use.names = FALSE)))
  stage_is_parallel <- vapply(topology, is.list, logical(1))

  G <- matrix(stats::rnorm(N * P), nrow = N)
  colnames(G) <- paste0("G", seq_len(P))
  CoG <- NULL
  if (n_CoG > 0) {
    CoG <- matrix(stats::rnorm(N * n_CoG), nrow = N)
    colnames(CoG) <- paste0("CoG", seq_len(n_CoG))
  }
  CoY <- NULL
  if (n_CoY > 0) {
    CoY <- matrix(stats::rnorm(N * n_CoY), nrow = N)
    colnames(CoY) <- paste0("CoY", seq_len(n_CoY))
  }
  G_full <- if (is.null(CoG)) G else cbind(G, CoG)

  X     <- vector("list", n_stage)   # X[[s]][[l]] cluster per layer
  Z     <- vector("list", n_stage)   # matrix, or list of matrices if parallel
  mu    <- vector("list", n_stage)
  sigma <- vector("list", n_stage)
  beta  <- NULL
  delta <- vector("list", max(0L, n_stage - 1L))

  for (s in seq_len(n_stage)) {
    Ks <- stage_K[[s]]
    n_layer_s <- length(Ks)

    # Design for this stage: exposures at stage 1, upstream clusters after.
    if (s == 1L) {
      design <- cbind(1, G_full)
    } else {
      prev_dummies <- do.call(cbind, lapply(seq_along(stage_K[[s - 1L]]), function(l) {
        stats::model.matrix(
          ~ factor(X[[s - 1L]][[l]], levels = seq_len(stage_K[[s - 1L]][l]))
        )[, -1, drop = FALSE]
      }))
      design <- cbind(1, prev_dummies)
    }

    Xs <- vector("list", n_layer_s)
    Zs <- vector("list", n_layer_s)
    mus <- vector("list", n_layer_s)
    sigmas <- vector("list", n_layer_s)
    coefs <- vector("list", n_layer_s)

    for (l in seq_len(n_layer_s)) {
      K_sl <- Ks[l]
      b <- matrix(0, nrow = K_sl, ncol = ncol(design))
      for (k in seq.int(2, K_sl)) {
        if (s == 1L) {
          b[k, 2] <- effect * (k - 1)                     # first exposure only
          if (n_CoG > 0) b[k, P + 1L + seq_len(n_CoG)] <- coef_CoG
        } else {
          b[k, -1] <- effect * (k - 1)                    # all upstream dummies
        }
      }
      coefs[[l]] <- b
      Xs[[l]] <- sim_rcat(sim_softmax(design %*% t(b)))

      mus[[l]]    <- sim_default_mu(K_sl, M, sep)
      sigmas[[l]] <- sim_default_sigma(K_sl, M)
      Zl <- matrix(0, nrow = N, ncol = M)
      for (k in seq_len(K_sl)) {
        idx <- which(Xs[[l]] == k)
        if (length(idx)) Zl[idx, ] <- sim_rmvnorm(length(idx), mus[[l]][k, ], sigmas[[l]][, , k])
      }
      colnames(Zl) <- paste0("S", s, "L", l, "_", seq_len(M))
      Zs[[l]] <- Zl
    }

    X[[s]] <- Xs
    mu[[s]] <- mus
    sigma[[s]] <- sigmas
    if (s == 1L) beta <- coefs else delta[[s - 1L]] <- coefs
    # A parallel stage keeps its list of layers; an early stage is one matrix.
    Z[[s]] <- if (stage_is_parallel[s]) Zs else Zs[[1L]]
  }

  # Outcome: driven by the last stage, summed across its layers.
  last_K <- stage_K[[n_stage]]
  gamma <- lapply(seq_along(last_K), function(l) seq_len(last_K[l]) - 1)
  eta_y <- rowSums(vapply(seq_along(last_K),
                          function(l) gamma[[l]][X[[n_stage]][[l]]], numeric(N)))
  if (n_CoY > 0) eta_y <- eta_y + as.numeric(CoY %*% rep(coef_CoY, n_CoY))
  Y <- if (family == "normal") stats::rnorm(N, eta_y, 1) else
       stats::rbinom(N, 1, stats::plogis(eta_y))
  Y <- as.matrix(Y); colnames(Y) <- "outcome"

  Z_complete <- Z

  # ---- missingness, applied per layer, preserving the nesting --------------
  if (missing != "none") {
    w <- if (mechanism == "MAR") {
      stats::plogis(0.9 * G[, 1] + 0.9 * as.numeric(scale(as.numeric(Y))))
    } else rep(1, N)
    blank <- function(Za) {
      if (missing %in% c("listwise", "mixed")) {
        rows <- sample.int(N, max(1L, round(miss_ratio * N)), prob = w / sum(w))
        Za[rows, ] <- NA
      }
      if (missing %in% c("sporadic", "mixed")) {
        pool <- which(!is.na(Za[, 1]))
        n_sp <- min(max(1L, round(miss_ratio * length(pool) * M * 0.5)),
                    length(pool) * (M - 1L))
        cells <- cbind(sample(pool, n_sp, replace = TRUE),
                       sample.int(M, n_sp, replace = TRUE))
        for (i in seq_len(nrow(cells))) {
          if (sum(is.na(Za[cells[i, 1], ])) < M - 1L) Za[cells[i, 1], cells[i, 2]] <- NA
        }
      }
      Za
    }
    Z <- lapply(Z, function(zs) if (is.list(zs)) lapply(zs, blank) else blank(zs))
  }

  list(
    G = G, Z = Z, Y = Y, CoG = CoG, CoY = CoY,
    Z_complete = Z_complete,
    K = topology,
    truth = list(
      model = "serial", family = family, topology = topology,
      stage_K = stage_K, stage_is_parallel = stage_is_parallel,
      N = N, P = P, M = M, n_stage = n_stage,
      beta = beta, delta = delta, mu = mu, sigma = sigma,
      gamma = gamma, X = X,
      n_CoG = n_CoG, n_CoY = n_CoY, coef_CoG = coef_CoG, coef_CoY = coef_CoY,
      mechanism = mechanism, missing = missing, miss_ratio = miss_ratio
    )
  )
}

#' Assemble the top-level \code{em_control} for a LUCID in serial fit
#'
#' A serial fit runs no EM loop of its own: it is a sequence of sub-model
#' fits, each with its own convergence record. Reporting only the stopping
#' controls here left a user with no top-level way to see that a stage had
#' exhausted \code{max_itr} -- they had to reach into
#' \code{submodel[[i]]$em_control} themselves. The settings are therefore
#' returned alongside diagnostics aggregated over the sub-models, with
#' \code{converged} \code{TRUE} only when every stage converged.
#'
#' @param submodel List of the fitted sub-models, one per stage.
#' @param tol,max_itr,max_tot.itr The stopping controls passed to every stage.
#' @return A list: \code{tol}, \code{max_itr}, \code{max_tot.itr} (the
#'   settings), \code{converged} (all stages converged), \code{n_iter} (total
#'   across stages), \code{submodel_converged}, \code{submodel_n_iter}, and
#'   \code{submodel_loglik_trace} (each stage's own log-likelihood trace, as
#'   a list rather than a single vector, since a serial fit has no EM loop of
#'   its own to trace).
#' @noRd
serial_em_control <- function(submodel, tol, max_itr, max_tot.itr) {
  sub_conv <- vapply(submodel,
                     function(s) isTRUE(s$em_control$converged),
                     logical(1))
  sub_iter <- vapply(submodel,
                     function(s) {
                       it <- s$em_control$n_iter
                       if (is.null(it)) NA_integer_ else as.integer(it)[1]
                     },
                     integer(1))
  # Each stage keeps its own loglik_trace; a serial fit has no EM loop of its
  # own to trace, so this is a list of each stage's trace rather than a
  # single vector the way early/parallel's em_control$loglik_trace is.
  sub_trace <- lapply(submodel, function(s) s$em_control$loglik_trace)
  nm <- names(submodel)
  if (!is.null(nm)) {
    names(sub_conv) <- nm
    names(sub_iter) <- nm
    names(sub_trace) <- nm
  }
  list(tol = tol,
       max_itr = max_itr,
       max_tot.itr = max_tot.itr,
       converged = all(sub_conv),
       n_iter = sum(sub_iter, na.rm = TRUE),
       submodel_converged = sub_conv,
       submodel_n_iter = sub_iter,
       submodel_loglik_trace = sub_trace)
}

#' @title Fit LUCID models with one or multiple omics layers
#' @description EM algorithm to estimate LUCID with one or multiple omics layers
#' @param lucid_model Specifying LUCID model, "early" for early integration, "parallel" for lucid in parallel,
#' "serial" for lucid in serial
#' @param G an N by P matrix representing exposures
#' @param Z Omics data, if "early", an N by M matrix; If "parallel", a list, each element i is a matrix with N rows and P_i features;
#' If "serial", a list, each element i is a matrix with N rows and p_i features or a list with two or more matrices with N rows and a certain number of features
#' @param Y a length N vector
#' @param CoG an N by V matrix representing covariates to be adjusted for G -> X
#' @param CoY an N by K matrix representing covariates to be adjusted for X -> Y
#' @param K Number of latent clusters. If "early", an integer greater or equal to 2; If "parallel", an integer vector, same length as Z, with each element being an integer greater or equal to 2;
#' If "serial", a list, each element is either an integer like that for "early" or an list of integers like that for "parallel", same length as Z
#' @param init_omic.data.model a vector of strings specifies the geometric model of omics
#' data. If NULL, See more in ?mclust::mclustModelNames
#' @param useY logical, if TRUE, EM algorithm fits a supervised LUCID; otherwise
#' unsupervised LUCID.
#' @param tol stopping criterion for the EM algorithm
#' @param max_itr Maximum iterations of the EM algorithm. If the EM algorithm iterates
#' more than max_itr without converging, the EM algorithm is forced to stop.
#' @param max_tot.itr Max number of total iterations for \code{estimate_lucid} function.
#' \code{estimate_lucid} may conduct EM algorithm for multiple times if the algorithm
#' fails to converge.
#' @param Rho_G A scalar. This parameter is the LASSO penalty to regularize
#' exposure coefficients in the G-to-X model. \code{CoG} adjustment covariates
#' are included unpenalized. If user wants to tune the penalty, use the wrapper
#' function \code{lucid}. Penalty tuning is supported for "early" and
#' "parallel". For "serial", only scalar penalty inputs are supported.
#' @param Rho_Z_Mu A scalar. This parameter is the LASSO penalty to
#' regularize cluster-specific means for omics data (Z). If user wants to tune the
#' penalty, use the wrapper function \code{lucid}. Penalty tuning is supported
#' for "early" and "parallel". For "serial", only scalar penalty inputs are
#' supported.
#' @param Rho_Z_Cov A scalar. This parameter is the graphical LASSO
#' penalty to estimate sparse cluster-specific variance-covariance matrices for omics
#' data (Z). If user wants to tune the penalty, use the wrapper function \code{lucid}.
#' Penalty tuning is supported for "early" and "parallel". For "serial", only
#' scalar penalty inputs are supported.
#' @param family The distribution of the outcome
#' @param seed Random seed to initialize the EM algorithm
#' @param init_impute Method to initialize the imputation of missing values in
#' LUCID. \code{lod} (the default) initializes the imputation via replacing
#' missing values by LOD / sqrt(2), where LOD is determined by the minimum of
#' each variable in omics data; \code{mix} uses \code{mclust::imputeData} to
#' implement EM Algorithm for Unrestricted General Location Model via the
#' \code{mix} package to impute the missing values in omics data. \code{mix}
#' is archived on CRAN and must be installed manually (e.g. from the CRAN
#' Archive) to use this option; a request for \code{init_impute = "mix"}
#' without \code{mix} installed raises an informative error.
#' @param init_par For "early", an interface to initialize EM algorithm, if mclust,
#' initiate the parameters using the \code{mclust} package, if random, initiate the parameters
#' by drawing from a uniform distribution;
#' For "parallel", mclust is the default for quick convergence;
#' For "serial", each sub-model follows the above depending on it is a "early" or "parallel"
#' @param verbose Logging level for fitting progress. If \code{FALSE}, concise
#' start/finish status lines are printed. If \code{TRUE}, detailed iteration-level
#' traces (including log-likelihood updates) are printed.
#' @param n_starts Number of independent random starts for the EM algorithm
#' (default 1). The EM algorithm converges only to a local optimum, so with
#' \code{n_starts > 1} the model is fitted from that many starting points and the
#' fit with the highest observed-data log-likelihood is returned. Per-start
#' log-likelihoods are recorded in \code{em_control$start_loglik}, which is worth
#' inspecting: a wide spread indicates the likelihood surface is multi-modal and
#' that a single start would have been unreliable.
#'
#' @import mclust
#' @import stats
#' @import utils
#' @import glasso
#' @import glmnet
#' @return An object of class \code{early_lucid}, \code{lucid_parallel} or
#' \code{lucid_serial} according to \code{lucid_model}. All three are lists;
#' the components common to every fit are:
#'
#' \describe{
#'   \item{res_Beta}{Estimates of the exposure-to-cluster (G -> X) association.
#'     For "early", a \code{K} by \code{(1 + P + V)} matrix of multinomial
#'     logistic coefficients, cluster 1 as reference. For "parallel" and
#'     "serial", a list holding the fitted object and the coefficient matrix
#'     per layer or stage.}
#'   \item{res_Mu}{Cluster-specific omics means (the mu of X -> Z). A
#'     \code{K} by \code{M} matrix for "early"; a list by layer or stage
#'     otherwise.}
#'   \item{res_Sigma}{Cluster-specific omics variance-covariance matrices (the
#'     sigma of X -> Z). A list of \code{K} matrices for "early"; a list by
#'     layer or stage otherwise.}
#'   \item{res_Gamma}{Estimates of the cluster-to-outcome (X -> Y) association,
#'     holding \code{beta} (absolute cluster levels), the reference-coded
#'     \code{cluster_effect} contrasts printed by \code{summary()}, any
#'     \code{covariate} coefficients, the residual \code{sigma} for a normal
#'     outcome, and the \code{parameterization} used.}
#'   \item{inclusion.p}{Posterior probability of cluster membership for each
#'     observation, \eqn{r_{ij}} of Eq 3. An \code{N} by \code{K} matrix for
#'     "early"; a list by layer or stage otherwise.}
#'   \item{K}{Number of latent clusters: an integer for "early", a list of
#'     integers for "parallel" and "serial".}
#'   \item{var.names}{Names of the G, Z and Y variables, as
#'     \code{list(Gnames, Znames, Ynames)}.}
#'   \item{init_omic.data.model}{The \code{mclust} geometric model used for the
#'     omics covariances.}
#'   \item{family}{Outcome distribution, "normal" or "binary".}
#'   \item{useY}{Whether the outcome was used in fitting, i.e. whether the model
#'     is supervised.}
#'   \item{Z}{The omics data. For "early" and "parallel" this is the data the
#'     model was fitted to, so sporadically missing cells hold their imputed
#'     values and listwise-missing rows remain \code{NA}. For "serial" it is the
#'     omics data \emph{as supplied}, still containing every missing value:
#'     imputation happens inside each stage, so the imputed omics for stage
#'     \code{i} are in \code{submodel[[i]]$Z}.}
#'   \item{init_impute}{The imputation method used to initialize missing omics
#'     values.}
#'   \item{init_par}{The parameter-initialization method used.}
#'   \item{Rho}{The penalties actually applied, as
#'     \code{list(Rho_G, Rho_Z_Mu, Rho_Z_Cov)}.}
#'   \item{missing_summary}{How much omics data was missing and in what pattern,
#'     using the taxonomy of the incomplete-omics extension: \code{complete_rows},
#'     \code{listwise_rows} (a whole omics layer missing for that subject) and
#'     \code{sporadic_rows} (some features missing), the corresponding
#'     proportions, and cell-level counts \code{total_missing_cells},
#'     \code{sporadic_missing_cells} and \code{prop_total_missing_cells}. A
#'     list by layer for "parallel"; for "serial", \code{n_stages} and a
#'     per-stage breakdown.}
#'   \item{em_control}{The stopping controls used (\code{tol}, \code{max_itr},
#'     \code{max_tot.itr}), which bootstrap refits reuse, together with
#'     convergence diagnostics. For "early" and "parallel" these are
#'     \code{converged} (\code{FALSE}, with a warning, if the fit exhausted
#'     \code{max_itr} without meeting \code{tol}), \code{n_iter},
#'     \code{n_restart}, \code{loglik_trace}, \code{n_starts} and
#'     \code{n_starts_ok}. A "serial" fit runs no EM loop of its own, so it
#'     reports \code{converged} (\code{TRUE} only if every sub-model converged),
#'     \code{n_iter} (the total across sub-models), the per-sub-model
#'     \code{submodel_converged} and \code{submodel_n_iter}, and
#'     \code{submodel_loglik_trace} (each sub-model's own \code{loglik_trace},
#'     by stage); the full diagnostics for a stage are in
#'     \code{submodel[[i]]$em_control}. For all three model types, the
#'     log-likelihood trace is checked for monotonicity as it is recorded: a
#'     decrease beyond a small, majorization-step-aware slack triggers a
#'     \code{warning()} naming the iteration (or stage) and the two values,
#'     since that indicates a numerical problem rather than expected EM
#'     behaviour.}
#' }
#'
#' The remaining components appear only for some model types:
#'
#' \describe{
#'   \item{likelihood}{The observed-data log-likelihood at the returned
#'     estimates, present for all three model types. For "early" and
#'     "parallel" this is the single joint log-likelihood from the EM fit. A
#'     "serial" fit runs no joint EM loop -- it is a sequence of conditionally
#'     fitted stages -- so its \code{likelihood} is instead the \emph{sum} of
#'     each stage's own log-likelihood (matching \code{cal_loglik_serial()});
#'     the individual per-stage values are in \code{submodel[[i]]$likelihood}.}
#'   \item{select}{Feature-selection indicators, present for all three model
#'     types. \code{select$selectG} and \code{select$selectZ} are logical
#'     vectors of retained exposures and omics features. For "parallel",
#'     \code{select$selectG} is the exposure-wise union across layers (selected
#'     in at least one layer), \code{select$selectG_layer} holds the per-layer
#'     exposure selection, and \code{select$selectZ} holds the per-layer omics
#'     selection (a list with at least one layer's selection reported). For
#'     "serial", \code{select} is \strong{stage 1's own selection only}
#'     (\code{submodel[[1]]$select}) -- stage 1 is the only stage whose G is
#'     the user's actual exposures; every later stage's G is the previous
#'     stage's posterior cluster-membership probabilities, so \code{Rho_G} is
#'     always 0 there and that stage's \code{selectG} is not a meaningful
#'     exposure-selection result (its \code{selectZ} still is). The complete
#'     per-stage selection record remains available at
#'     \code{submodel[[i]]$select} for every \code{i}.}
#'   \item{N}{Number of observations. Present for "parallel" and "serial"; for
#'     "early", use \code{nrow(fit$inclusion.p)}.}
#'   \item{z}{The E-step responsibilities over the joint cluster configuration
#'     across layers, before they are marginalised into \code{inclusion.p}.
#'     Present for "parallel" only -- nothing in the package reads this field,
#'     and there is no plan to extend it to "serial".}
#'   \item{res_Delta}{Estimates of the between-stage cluster transition
#'     associations, one element per transition (length \code{n_stages - 1}).
#'     Element \code{i} is the coefficient object of stage \code{i + 1} fitted
#'     with stage \code{i}'s cluster assignment in place of the exposures, so it
#'     has the same structure as \code{res_Beta} for that stage. Present for
#'     "serial" only.}
#'   \item{submodel}{The fitted sub-models, one per stage, each itself an
#'     \code{early_lucid} or \code{lucid_parallel} object. Present for "serial"
#'     only.}
#' }
#' @examples
#' i <- 1008
#' set.seed(i)
#' G <- matrix(rnorm(500), nrow = 100)
#' Z1 <- matrix(rnorm(1000), nrow = 100)
#' Z2 <- matrix(rnorm(1000), nrow = 100)
#' Z3 <- matrix(rnorm(1000), nrow = 100)
#' Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3)
#' Y <- rnorm(100)
#' CoY <- matrix(rnorm(200), nrow = 100)
#' CoG <- matrix(rnorm(200), nrow = 100)
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2, 2, 2),
#' lucid_model = "serial",
#' family = "normal",
#' seed = i,
#' CoG = CoG, CoY = CoY,
#' useY = TRUE,
#' max_itr = 20, max_tot.itr = 50)
#' @export
#'
estimate_lucid <- function(lucid_model = c("early", "parallel","serial"),
                      G, Z, Y, CoG = NULL, CoY = NULL, K,
                      init_omic.data.model = "EEV",
                      useY = TRUE,
                      tol = 1e-3,
                      max_itr = 1e3,
                      max_tot.itr = 1e4,
                      Rho_G = 0,
                      Rho_Z_Mu = 0,
                      Rho_Z_Cov = 0,
                      family = c("normal", "binary"),
                      seed = 123,
                      init_impute = c("lod", "mix"),
                      init_par = c("mclust", "random"),
                      verbose = FALSE,
                      n_starts = 1L) {
  family <- normalize_family_label(family)
  # Resolve here rather than leaving it to est_lucid(). The serial path assembles
  # its own result list and stored whatever was passed, so a caller who took the
  # default got the unevaluated c("mix", "lod") / c("mclust", "random") back on
  # the fit -- a record of the choice that was never made, and one that a
  # bootstrap refit reading these fields would carry forward.
  init_impute <- match.arg(init_impute)
  init_par <- match.arg(init_par)
  # Resolve lucid_model here too, for the same reason: a caller who omitted it
  # (relying on the documented "early" default) previously had the unevaluated
  # 3-element default vector forwarded to est_lucid() as-is, which has its own,
  # different default and no way to tell "not supplied" from "supplied the
  # wrong length" -- match.arg() inside est_lucid() then rejected it outright,
  # so estimate_lucid() could not actually be called without naming
  # lucid_model explicitly, despite its signature advertising a default.
  lucid_model <- match.arg(lucid_model)
  if (lucid_model == "early" || lucid_model == "parallel"){
    # ========================== Early Integration ==========================
    # ========================== LUCID IN PARALLEL ==========================
    results <- est_lucid(lucid_model = lucid_model,
                         G = G,
                         Z = Z,
                         Y = Y,
                         CoG = CoG,
                         CoY = CoY, K = K,
                         init_omic.data.model = init_omic.data.model,
                         useY = useY,
                         tol = tol,
                         max_itr = max_itr,
                         max_tot.itr = max_tot.itr,
                         Rho_G = Rho_G,
                         Rho_Z_Mu = Rho_Z_Mu,
                         Rho_Z_Cov = Rho_Z_Cov,
                         family = family,
                         seed = seed,
                         init_impute = init_impute,
                         init_par = init_par,
                         verbose = verbose,
                         n_starts = n_starts)
    return(results)
  }else{
    # ========================== LUCID IN Serial ==========================
    if (!is.list(K)) {
      if (is.numeric(K) && length(K) > 0) {
        # Backward-compatible support: serial with only early stages may pass K as numeric vector.
        K <- as.list(as.numeric(K))
      } else {
        stop("For LUCID in Serial, K should be a non-empty list!")
      }
    }
    if (length(K) == 0) {
      stop("For LUCID in Serial, K should be a non-empty list!")
    }
    if (is.null(Z)) {
      stop("Input data 'Z' is missing")
    }
    if (!is.list(Z)) {
      stop("Input data 'Z' should be a list for LUCID in Serial!")
    }
    if (length(Z) != length(K)) {
      stop("Z and K should be two lists of the same length for LUCID in Serial!")
    }

    # Validate G/Y/CoG/CoY completeness here, at the top level, rather than
    # letting a missing value surface only once stage 1's nested est_lucid()
    # re-runs this same check -- otherwise the error reports as coming from
    # inside a sub-fit instead of from the call the user actually made.
    check_complete_input(G, "G")
    check_complete_input(Y, "Y")
    check_complete_input(CoG, "CoG")
    check_complete_input(CoY, "CoY")

    normalize_serial_block <- function(z_block, k_block, block_id = "root") {
      if (is.list(k_block)) {
        if (!is.list(z_block)) {
          stop(paste0("For LUCID in Serial, input data 'Z' must match 'K' structure. Nested K at ",
                      block_id, " requires nested Z list."))
        }
        if (length(z_block) != length(k_block)) {
          stop(paste0("For LUCID in Serial, nested Z length must equal nested K length at ",
                      block_id, "."))
        }
        out <- vector("list", length(z_block))
        for (j in seq_along(z_block)) {
          out[[j]] <- normalize_serial_block(z_block[[j]], k_block[[j]],
                                             block_id = paste0(block_id, ".", j))
        }
        names(out) <- names(z_block)
        return(out)
      }

      if (!is.numeric(k_block) || length(k_block) != 1 || is.na(k_block) || k_block < 2) {
        stop(paste0("For estimate_lucid(..., lucid_model = 'serial'), each non-list K entry ",
                    "must be a single integer >= 2. Invalid entry at ", block_id, "."))
      }

      z_mat <- as.matrix(z_block)
      if (!is.numeric(z_mat)) {
        stop(paste0("For LUCID in Serial, Z block at ", block_id, " must be numeric."))
      }
      z_mat
    }

    for (i in seq_along(K)) {
      Z[[i]] <- normalize_serial_block(Z[[i]], K[[i]], block_id = paste0("stage", i))
    }

    extract_nonref_pip <- function(model_obj, n_obs) {
      if (inherits(model_obj, "early_lucid")) {
        p <- as.matrix(model_obj$inclusion.p)
        if (ncol(p) <= 1) {
          return(matrix(numeric(0), nrow = n_obs))
        }
        return(p[, -1, drop = FALSE])
      }
      if (inherits(model_obj, "lucid_parallel")) {
        p_list <- model_obj$inclusion.p
        p_nonref <- lapply(p_list, function(p) {
          p <- as.matrix(p)
          if (ncol(p) <= 1) {
            matrix(numeric(0), nrow = nrow(p))
          } else {
            p[, -1, drop = FALSE]
          }
        })
        if (length(p_nonref) == 0) {
          return(matrix(numeric(0), nrow = n_obs))
        }
        return(do.call(cbind, p_nonref))
      }
      stop("Unsupported submodel class in serial pipeline.")
    }

    fit_serial_stage <- function(stage_idx,
                                 stage_model,
                                 G_stage,
                                 Z_stage,
                                 Y_stage,
                                 CoG_stage,
                                 CoY_stage,
                                 K_stage,
                                 useY_stage,
                                 family_stage) {
      G_stage <- as.matrix(G_stage)
      if (!is.numeric(G_stage)) {
        stop("Serial stage input G must be numeric.")
      }

      # Stage 1's G is the user's real exposures; every later stage's G is the
      # *previous* stage's posterior cluster-membership probabilities
      # (stage_G <- post.p below), not an exposure. Penalizing/selecting on
      # that pseudo-G has no interpretation as exposure selection, so Rho_G is
      # forced to 0 there unconditionally -- not only when it happens to have
      # fewer than 2 columns.
      rho_g_stage <- Rho_G
      if (rho_g_stage > 0 && stage_idx > 1) {
        rho_g_stage <- 0
        if (verbose) {
          cat(sprintf("Stage %d: Rho_G forced to 0 -- this stage's G is the previous stage's cluster assignment, not an exposure.\n",
                      stage_idx))
        }
      }

      if (isTRUE(verbose)) {
        return(est_lucid(
          lucid_model = stage_model,
          G = G_stage,
          Z = Z_stage,
          Y = Y_stage,
          CoG = CoG_stage,
          CoY = CoY_stage,
          K = K_stage,
          init_omic.data.model = init_omic.data.model,
          useY = useY_stage,
          tol = tol,
          max_itr = max_itr,
          max_tot.itr = max_tot.itr,
          Rho_G = rho_g_stage,
          Rho_Z_Mu = Rho_Z_Mu,
          Rho_Z_Cov = Rho_Z_Cov,
          family = family_stage,
          seed = seed + stage_idx * 1900,
          init_impute = init_impute,
          init_par = init_par,
          verbose = TRUE,
          relabel_by_outcome = stage_idx == n_stage,
          n_starts = n_starts
        ))
      }

      stage_fit <- NULL
      invisible(capture.output(
        stage_fit <- est_lucid(
          lucid_model = stage_model,
          G = G_stage,
          Z = Z_stage,
          Y = Y_stage,
          CoG = CoG_stage,
          CoY = CoY_stage,
          K = K_stage,
          init_omic.data.model = init_omic.data.model,
          useY = useY_stage,
          tol = tol,
          max_itr = max_itr,
          max_tot.itr = max_tot.itr,
          Rho_G = rho_g_stage,
          Rho_Z_Mu = Rho_Z_Mu,
          Rho_Z_Cov = Rho_Z_Cov,
          family = family_stage,
          seed = seed + stage_idx * 1900,
          init_impute = init_impute,
          init_par = init_par,
          verbose = FALSE,
          relabel_by_outcome = stage_idx == n_stage,
          n_starts = n_starts
        )
      ))
      stage_fit
    }

    n_stage <- length(K)
    post.p.list <- vector(mode = "list", length = n_stage)
    res.mu.list <- vector(mode = "list", length = n_stage)
    res.sigma.list <- vector(mode = "list", length = n_stage)
    res.delta.list <- vector(mode = "list", length = max(0, n_stage - 1))
    Znames <- vector(mode = "list", length = n_stage)
    submodel <- vector(mode = "list", length = n_stage)
    missing_by_stage <- vector(mode = "list", length = n_stage)
    has_penalty <- function(model_obj) {
      if (is.null(model_obj$Rho)) return(FALSE)
      any(unlist(model_obj$Rho[c("Rho_G", "Rho_Z_Mu", "Rho_Z_Cov")]) != 0)
    }
    stage_selection_msg <- function(model_obj) {
      if (inherits(model_obj, "early_lucid")) {
        g_sel <- sum(model_obj$select$selectG)
        g_tot <- length(model_obj$select$selectG)
        z_sel <- sum(model_obj$select$selectZ)
        z_tot <- length(model_obj$select$selectZ)
        return(sprintf("Selected G: %d/%d; Selected Z: %d/%d.", g_sel, g_tot, z_sel, z_tot))
      }
      if (inherits(model_obj, "lucid_parallel")) {
        g_sel <- sum(model_obj$select$selectG)
        g_tot <- length(model_obj$select$selectG)
        z_sel <- vapply(model_obj$select$selectZ, function(x) {
          if (is.null(dim(x))) sum(x) else sum(colSums(x) > 0)
        }, numeric(1))
        z_tot <- vapply(model_obj$select$selectZ, function(x) {
          if (is.null(dim(x))) length(x) else ncol(x)
        }, numeric(1))
        return(sprintf("Selected G: %d/%d; Selected Z by layer: %s.",
                       g_sel, g_tot, paste0(z_sel, "/", z_tot, collapse = ", ")))
      }
      ""
    }

    if (!isTRUE(verbose)) {
      cat(sprintf("Fitting LUCID serial model (%d stages)...\n", n_stage))
    }

    post.p <- NULL
    for (stage_idx in seq_len(n_stage)) {
      if (verbose) {
        cat(sprintf("Fitting LUCID serial model (Stage %d/%d)...\n", stage_idx, n_stage))
      }

      is_last <- (stage_idx == n_stage)
      stage_model <- if (is.list(K[[stage_idx]])) "parallel" else "early"
      stage_K <- if (is.list(K[[stage_idx]])) {
        as.numeric(unlist(K[[stage_idx]], use.names = FALSE))
      } else {
        as.numeric(K[[stage_idx]])
      }
      stage_Y <- if (is_last) Y else runif(nrow(G))
      stage_family <- if (is_last) family else "normal"
      stage_useY <- if (is_last) useY else FALSE
      stage_CoY <- if (is_last) CoY else NULL
      stage_CoG <- if (stage_idx == 1) CoG else NULL
      stage_G <- if (stage_idx == 1) G else post.p

      temp_model <- fit_serial_stage(
        stage_idx = stage_idx,
        stage_model = stage_model,
        G_stage = stage_G,
        Z_stage = Z[[stage_idx]],
        Y_stage = stage_Y,
        CoG_stage = stage_CoG,
        CoY_stage = stage_CoY,
        K_stage = stage_K,
        useY_stage = stage_useY,
        family_stage = stage_family
      )
      if (!isTRUE(verbose)) {
        if (has_penalty(temp_model)) {
          cat(sprintf("  Stage %d/%d (%s) finished: log-likelihood = %.3f. %s\n",
                      stage_idx, n_stage, stage_model, temp_model$likelihood,
                      stage_selection_msg(temp_model)))
        } else {
          cat(sprintf("  Stage %d/%d (%s) finished: log-likelihood = %.3f.\n",
                      stage_idx, n_stage, stage_model, temp_model$likelihood))
        }
      }

      post.p.list[[stage_idx]] <- temp_model$inclusion.p
      res.mu.list[[stage_idx]] <- temp_model$res_Mu
      res.sigma.list[[stage_idx]] <- temp_model$res_Sigma
      Znames[[stage_idx]] <- temp_model$var.names$Znames
      submodel[[stage_idx]] <- temp_model
      missing_by_stage[[stage_idx]] <- temp_model$missing_summary

      if (stage_idx == 1) {
        res_Beta <- temp_model$res_Beta
        Gnames <- temp_model$var.names$Gnames
      } else {
        res.delta.list[[stage_idx - 1]] <- temp_model$res_Beta
      }

      if (is_last) {
        res_Gamma <- temp_model$res_Gamma
        Ynames <- temp_model$var.names$Ynames
      } else {
        post.p <- extract_nonref_pip(temp_model, n_obs = nrow(G))
        if (ncol(post.p) == 0) {
          stop(paste0("Sub Model ", stage_idx,
                      " produced no non-reference cluster probabilities to pass forward."))
        }
      }
    }

    if(verbose){
      cat("Success: LUCID serial model constructed!", "\n\n")
    } else {
      cat("Finished LUCID serial model.\n")
    }
    serial_missing_summary <- list(
      n_stages = n_stage,
      stage = missing_by_stage
    )
    results <- list(res_Beta = res_Beta,
                    res_Mu = res.mu.list,
                    res_Sigma = res.sigma.list,
                    res_Delta = res.delta.list,
                    res_Gamma = res_Gamma,
                    K = K,
                    N = nrow(G),
                    var.names =list(Gnames = Gnames,
                                    Znames = Znames,
                                    Ynames = Ynames),
                    init_omic.data.model =  init_omic.data.model,
                    inclusion.p = post.p.list,
                    family = family,
                    useY = useY,
                    Z = Z,
                    init_impute = init_impute,
                    init_par = init_par,
                    submodel = submodel,
                    missing_summary = serial_missing_summary,
                    Rho = list(Rho_G = Rho_G,
                               Rho_Z_Mu = Rho_Z_Mu,
                               Rho_Z_Cov = Rho_Z_Cov),
                    em_control = serial_em_control(submodel, tol = tol,
                                                   max_itr = max_itr,
                                                   max_tot.itr = max_tot.itr)
    )
    class(results) <- c("lucid_serial")

    # `likelihood`: a serial chain is a sequence of conditional fits, not one
    # joint EM, so there is no single joint log-likelihood to report. We
    # instead surface the sum of each stage's own likelihood, matching
    # cal_loglik_serial() (R/summary.R) exactly -- this makes `likelihood` a
    # field present on all three model types.
    results$likelihood <- cal_loglik_serial(results)

    # `select`: only stage 1's G is the user's real exposures -- every later
    # stage's G is the previous stage's cluster-membership probabilities
    # (Rho_G is forced to 0 there; see fit_serial_stage() above), so only
    # stage 1's selectG means "which exposures were selected." This mirrors
    # how res_Beta above is already stage-1-scoped. The complete per-stage
    # record (including each stage's selectZ, which is meaningful throughout)
    # remains at submodel[[i]]$select for every i, as consumed today by
    # has_unselected_feature_serial() (R/boot_lucid.R).
    results$select <- submodel[[1]]$select

    # No serial `z`: `z` (the joint E-step responsibility array, present on
    # parallel fits only) has no downstream reader anywhere in the package,
    # and the only code that ever read it -- reorder_lucid()/reorder_z(),
    # themselves already dead and removed -- is gone. There is nothing this
    # field would serve even if resurrected.

    return(results)
  }
}


# =============================================================================
# est_lucid(): the EM workhorse for "early" and "parallel" (merged in from the
# former em.R). estimate_lucid() above dispatches to this for early/parallel and
# drives its own stage loop for serial; lucid.R and boot_lucid.R also call this
# directly (bypassing estimate_lucid()), so it is not exclusively the wrapper
# above's private engine, just its closest and most frequent caller.
# =============================================================================

#' EM workhorse for the "early" and "parallel" LUCID models
#'
#' Runs the EM algorithm to convergence for a single early- or
#' parallel-integration fit. This is the actual estimation engine:
#' \code{estimate_lucid()} dispatches to it directly for early/parallel and
#' calls it once per stage for serial; \code{lucid()} and \code{boot_lucid()}
#' also call it directly, bypassing \code{estimate_lucid()} entirely.
#'
#' \code{n_starts > 1} runs \code{n_starts} independent EM initializations
#' (via a recursive self-call with \code{n_starts = 1}) and keeps the fit
#' with the highest observed-data log-likelihood, since EM only guarantees a
#' local optimum and a single start can silently return a poor one.
#'
#' @param lucid_model "early" or "parallel".
#' @param G,Z,Y,CoG,CoY Exposure, omics, outcome, and optional covariate data.
#' @param K Number of latent clusters (early: scalar; parallel: vector, one
#'   per layer).
#' @param init_omic.data.model The \code{mclust} geometric model for the
#'   omics covariances.
#' @param useY Whether to use the outcome in fitting (supervised vs.
#'   unsupervised).
#' @param tol,max_itr,max_tot.itr Convergence tolerance and iteration limits.
#' @param Rho_G,Rho_Z_Mu,Rho_Z_Cov LASSO penalties for exposure and omics
#'   selection.
#' @param family "normal" or "binary" outcome.
#' @param seed Random seed for initialization.
#' @param init_impute,init_par Initialization strategy for missing omics and
#'   for starting parameter values.
#' @param verbose Whether to print per-iteration progress.
#' @param relabel_by_outcome Whether to relabel clusters by outcome risk at
#'   the end of fitting (set \code{FALSE} for an upstream serial stage, which
#'   has no real outcome to relabel by).
#' @param n_starts Number of independent EM initializations to run, keeping
#'   the best.
#' @return An object of class \code{early_lucid} or \code{lucid_parallel} --
#'   see \code{estimate_lucid()}'s \code{@return} for the field-by-field
#'   description, since that is this function's own output, surfaced as-is.
#' @noRd
est_lucid <-
function (lucid_model = c("early", "parallel"), G, Z, Y, CoG = NULL,
    CoY = NULL, K, init_omic.data.model = "EEV", useY = TRUE,
    tol = 0.001, max_itr = 1000, max_tot.itr = 10000, Rho_G = 0,
    Rho_Z_Mu = 0, Rho_Z_Cov = 0, family = c("normal", "binary"),
    seed = 123, init_impute = c("lod", "mix"), init_par = c("mclust",
        "random"), verbose = FALSE, relabel_by_outcome = TRUE,
    n_starts = 1L)
{
    # ---- Multi-start selection (D10) -------------------------------------
    # Run `n_starts` independent EM initializations and keep the fit with the
    # highest observed-data log-likelihood.  EM is only guaranteed to reach a
    # local optimum, so a single start can silently return a poor solution.
    n_starts <- as.integer(n_starts)
    if (is.na(n_starts) || n_starts < 1L) {
        stop("'n_starts' should be a positive integer")
    }
    if (n_starts > 1L) {
        cl <- match.call()
        cl$n_starts <- 1L
        best <- NULL
        best_loglik <- -Inf
        n_ok <- 0L
        start_loglik <- rep(NA_real_, n_starts)
        for (s in seq_len(n_starts)) {
            cl$seed <- seed + (s - 1L) * 997L
            fit_s <- try(eval(cl, parent.frame()), silent = TRUE)
            if (inherits(fit_s, "try-error")) next
            n_ok <- n_ok + 1L
            start_loglik[s] <- fit_s$likelihood
            if (is.finite(fit_s$likelihood) && fit_s$likelihood > best_loglik) {
                best_loglik <- fit_s$likelihood
                best <- fit_s
            }
        }
        if (is.null(best)) {
            stop("All ", n_starts, " EM starts failed; try different starting seeds or 'init_par'")
        }
        best$em_control$n_starts <- n_starts
        best$em_control$n_starts_ok <- n_ok
        best$em_control$start_loglik <- start_loglik
        return(best)
    }
    init_impute <- match.arg(init_impute)
    init_par <- match.arg(init_par)
    family <- normalize_family_label(family)
    Select_G <- FALSE
    Select_Z <- FALSE
    if (Rho_G != 0) {
        Select_G <- TRUE
    }
    if (Rho_Z_Mu != 0 | Rho_Z_Cov != 0) {
        Select_Z <- TRUE
    }
    penalty_requested <- isTRUE(Select_G) || isTRUE(Select_Z)
    if (is.null(G)) {
        stop("Input data 'G' is missing")
    }
    else {
        if (!is.matrix(G)) {
            G <- as.matrix(G)
            if (!is.numeric(G)) {
                stop("Input data 'G' should be numeric; categorical variables should be transformed into dummies")
            }
        }
    }
    if (is.null(colnames(G))) {
        Gnames <- paste0("G", 1:ncol(G))
    }
    else {
        Gnames <- colnames(G)
    }
    colnames(G) <- Gnames
    if (is.null(Y)) {
        stop("Input data 'Y' is missing")
    }
    else {
        if (!is.matrix(Y)) {
            Y <- as.matrix(Y)
            if (!is.numeric(Y)) {
                stop("Input data 'Y' should be numeric; binary outcome should be transformed them into dummies")
            }
            if (ncol(Y) > 1) {
                stop("Only continuous 'Y' or binary 'Y' is accepted")
            }
        }
    }
    if (is.null(colnames(Y))) {
        Ynames <- "outcome"
    }
    else {
        Ynames <- colnames(Y)
    }
    colnames(Y) <- Ynames
    check_complete_input(G, "G")
    check_complete_input(Y, "Y")
    check_complete_input(CoG, "CoG")
    check_complete_input(CoY, "CoY")
    if (is_binary_family(family)) {
        if (!(all(Y %in% c(0, 1)))) {
            stop("Binary outcome should be coded as 0 and 1")
        }
    }
    CoGnames <- NULL
    if (!is.null(CoG)) {
        if (!is.matrix(CoG)) {
            CoG <- as.matrix(CoG)
            if (!is.numeric(CoG)) {
                stop("Input data 'CoG' should be numeric; categroical variables should be transformed into dummies")
            }
        }
        if (is.null(colnames(CoG))) {
            CoGnames <- paste0("CoG", 1:ncol(CoG))
        }
        else {
            CoGnames <- colnames(CoG)
        }
        colnames(CoG) <- CoGnames
    }
    CoYnames <- NULL
    if (!is.null(CoY)) {
        if (!is.matrix(CoY)) {
            CoY <- as.matrix(CoY)
            if (!is.numeric(CoY)) {
                stop("Input data 'CoY' should be numeric; categorical variables should be transformed into dummies")
            }
        }
        if (is.null(colnames(CoY))) {
            CoYnames <- paste0("CoY", 1:ncol(CoY))
        }
        else {
            CoYnames <- colnames(CoY)
        }
        colnames(CoY) <- CoYnames
    }
    # K >= 2 is already enforced when fitting through lucid()/tune_lucid(),
    # but est_lucid() is documented and callable directly, and previously
    # skipped this check itself -- a K = 1 call reached deep EM machinery
    # before failing with an obscure error instead of this clear one.
    check_K(K)
    if (match.arg(lucid_model) == "early") {
        family <- normalize_family_label(family)
        init_omic.data.model = init_omic.data.model
        if (is.null(Z)) {
            stop("Input data 'Z' is missing")
        }
        else {
            if (!is.matrix(Z)) {
                Z <- as.matrix(Z)
                if (!is.numeric(Z)) {
                  stop("Input data 'Z' should be numeric")
                }
            }
        }
        if (is.null(colnames(Z))) {
            Znames <- paste0("Z", 1:ncol(Z))
        }
        else {
            Znames <- colnames(Z)
        }
        N <- nrow(Y)
        dimG <- ncol(G)
        dimZ <- ncol(Z)
        dimCoG <- ifelse(is.null(CoG), 0, ncol(CoG))
        dimCoY <- ifelse(is.null(CoY), 0, ncol(CoY))
        G <- cbind(G, CoG)
        Gnames <- c(Gnames, CoGnames)
        family.list <- switch(family, normal = normal(K = K, 
            dimCoY), binary = binary(K = K, dimCoY))
        Mstep_Y <- family.list$f.maxY
        switch_Y <- family.list$f.switch
        na_pattern <- check_na(Z)
        missing_summary <- summarize_missing_stats(na_pattern, lucid_model = "early")
        if (!isTRUE(verbose)) {
            cat(sprintf("Fitting LUCID early model (K = %s)...\n", paste(K, collapse = ",")))
        }
        if (na_pattern$impute_flag) {
            if (init_impute == "mix") {
                if (!requireNamespace("mix", quietly = TRUE)) {
                  stop("init_impute = \"mix\" requires the 'mix' package, ",
                    "which is archived on CRAN; install it manually (e.g. ",
                    "from the CRAN Archive) or use init_impute = \"lod\" ",
                    "instead.", call. = FALSE)
                }
                if (verbose) {
                  cat("Intializing imputation of missing values in 'Z' via the mix package \n\n")
                }
                invisible(capture.output(Z <- mclust::imputeData(Z,
                  seed = seed)))
                Z[na_pattern$indicator_na == 3, ] <- NA
            }
            if (init_impute == "lod") {
                if (verbose) {
                  cat("Intializing imputation of missing values in 'Z' via LOD / sqrt(2) \n\n")
                }
                Z <- apply(Z, 2, fill_data_lod)
                colnames(Z) <- Znames
            }
            # D7: rows with no observed omics value must stay NA regardless of
            # initializer.  The LOD branch fills every missing cell, which left
            # fabricated values in fit$Z for participants who were never measured.
            Z[na_pattern$indicator_na == 3, ] <- NA
        }
        tot.itr <- 0
        convergence <- FALSE
        em_converged <- FALSE
        em_iter <- 0L
        loglik_trace <- numeric(0)
        n_restart <- 0L
        while (!convergence && tot.itr <= max_tot.itr) {
            if (tot.itr > 0) {
                seed <- seed + 10
                n_restart <- n_restart + 1L
            }
            set.seed(seed)
            loglik_trace <- numeric(0)
            res.beta <- matrix(data = runif(K * (dimG + dimCoG + 
                1)), nrow = K)
            res.beta[1, ] <- 0
            if (init_par == "mclust") {
                if (verbose) {
                  cat("Initialize LUCID with mclust based on inclusion probabilities given by mclust \n")
                }
                invisible(capture.output(mclust.fit <- Mclust(Z[na_pattern$indicator_na != 
                  3, ], G = K, modelNames = init_omic.data.model)))
                if (is.null(mclust.fit)) {
                  stop("mclust failed for specified model - please set init_omic.data.model to `NULL` to conduct automatic model selection ")
                }
                if (is.null(init_omic.data.model)) {
                  model.best <- mclust.fit$modelName
                }
                else {
                  model.best <- init_omic.data.model
                }
                res.mu <- t(mclust.fit$parameters$mean)
                res.sigma <- mclust.fit$parameters$variance$sigma
            }
            else {
                if (verbose) {
                  cat("Initialize LUCID with random values from uniform distribution \n")
                }
                if (is.null(init_omic.data.model)) {
                  model.best <- "EEV"
                  if (verbose) {
                    cat("GMM model for LUCID is not specified, 'EEV' model is used by default \n")
                  }
                }
                else {
                  model.best <- init_omic.data.model
                }
                res.mu <- matrix(runif(dimZ * K, min = -0.5, 
                  max = 0.5), nrow = K)
                res.sigma <- gen_cov_matrices(dimZ = dimZ, K = K)
            }
            res.gamma <- family.list$initial.gamma(K, dimCoY)
            if (verbose) {
                cat("Fitting Early Integration LUCID model", 
                  paste0("(", "K = ", K, ", Rho_G = ", Rho_G, 
                    ", Rho_Z_Mu = ", Rho_Z_Mu, ", Rho_Z_Cov = ", 
                    Rho_Z_Cov, ")"), "\n")
            }
            res.loglik <- -Inf
            itr <- 0
            restart_needed <- FALSE
            while (!convergence && itr <= max_itr) {
                itr <- itr + 1
                tot.itr <- tot.itr + 1
                check.gamma <- TRUE
                new.likelihood <- Estep_early(beta = res.beta, 
                  mu = res.mu, sigma = res.sigma, gamma = res.gamma, 
                  G = G, Z = Z, Y = Y, CoY = CoY, N = N, K = K, 
                  family.list = family.list, itr = itr, useY = useY, 
                  dimCoY = dimCoY, ind.na = na_pattern$indicator_na)
                res.r <- t(apply(new.likelihood, 1, lse_vec))
                if (!all(is.finite(res.r))) {
                  if (verbose) {
                    cat("iteration", itr, ": EM algorithm collapsed: invalid estiamtes due to over/underflow, try LUCID with another seed \n")
                  }
                  restart_needed <- TRUE
                  break
                }
                else {
                  if (isTRUE(verbose)) {
                    cat("iteration", itr, ": E-step finished.\n")
                  }
                }
                invisible(capture.output(new.beta <- Mstep_G(G = G, 
                  r = res.r, selectG = Select_G, penalty = Rho_G, 
                  dimG = dimG, dimCoG = dimCoG, K = K)))
                new.mu.sigma <- Mstep_Z(Z = Z, r = res.r, selectZ = Select_Z, 
                  penalty.mu = Rho_Z_Mu, penalty.cov = Rho_Z_Cov, 
                  model.name = model.best, K = K, ind.na = na_pattern$indicator_na, 
                  mu = res.mu)
                if (is.null(new.mu.sigma$mu)) {
                  if (verbose) {
                    cat("variable selection failed, try LUCID with another seed \n")
                  }
                  restart_needed <- TRUE
                  break
                }
                if (useY) {
                  new.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, 
                    K = K, CoYnames)
                  check.gamma <- all(is.finite(c(new.gamma$beta, new.gamma$sigma)))
                }
                if (na_pattern$impute_flag) {
                  # The I-step maximizes the weighted Gaussian density with the
                  # responsibilities held fixed, so those weights must be the
                  # posterior at the SAME parameters used for the density.
                  # Algorithm 1 reuses r from the E-step at Theta^(t) while the
                  # density already uses Theta^(t+1); that mismatch breaks the
                  # ascent guarantee and made the observed log-likelihood dip by
                  # small amounts late in the run.  Recompute r at Theta^(t+1)
                  # first, which restores monotonicity and leaves the converged
                  # solution unchanged (the two coincide at a fixed point).
                  istep.r <- t(apply(Estep_early(
                    beta = new.beta, mu = new.mu.sigma$mu,
                    sigma = new.mu.sigma$sigma,
                    gamma = if (useY) new.gamma else res.gamma,
                    G = G, Z = Z, Y = Y, CoY = CoY, N = N, K = K,
                    family.list = family.list, itr = itr, useY = useY,
                    dimCoY = dimCoY, ind.na = na_pattern$indicator_na
                  ), 1, lse_vec))
                  if (all(is.finite(istep.r))) {
                    Z <- Istep_Z(Z = Z, p = istep.r, mu = new.mu.sigma$mu,
                      sigma = new.mu.sigma$sigma, index = na_pattern$index,
                      lucid_model = "early")
                  } else {
                    Z <- Istep_Z(Z = Z, p = res.r, mu = new.mu.sigma$mu,
                      sigma = new.mu.sigma$sigma, index = na_pattern$index,
                      lucid_model = "early")
                  }
                }
                check.value <- all(is.finite(new.beta), is.finite(unlist(new.mu.sigma)), 
                  check.gamma)
                if (!check.value) {
                  if (verbose) {
                    cat("iteration", itr, ": Invalid estimates, try LUCID with another seed \n")
                  }
                  restart_needed <- TRUE
                  break
                }
                else {
                  res.beta <- new.beta
                  res.mu <- new.mu.sigma$mu
                  res.sigma <- new.mu.sigma$sigma
                  if (useY) {
                    res.gamma <- new.gamma
                  }
                  updated.likelihood <- Estep_early(
                    beta = res.beta, mu = res.mu, sigma = res.sigma,
                    gamma = res.gamma, G = G, Z = Z, Y = Y, CoY = CoY,
                    N = N, K = K, family.list = family.list, itr = itr + 1L,
                    useY = useY, dimCoY = dimCoY,
                    ind.na = na_pattern$indicator_na
                  )
                  res.r <- t(apply(updated.likelihood, 1, lse_vec))
                  raw.loglik <- observed_loglik(updated.likelihood)
                  # `new.loglik` is the PENALIZED objective -- what the penalized
                  # EM actually ascends, so it drives the convergence/dip checks
                  # below. `loglik_trace` and the stored `$likelihood` are the
                  # unpenalized observed-data log-likelihood (`raw.loglik`), so
                  # the trace has one meaning across model types and its last
                  # value equals `$likelihood`.
                  new.loglik <- raw.loglik - early_penalty_value(
                    res.beta, res.mu, res.sigma, Rho_G, Rho_Z_Mu,
                    Rho_Z_Cov, dimG
                  )
                  if (isTRUE(verbose)) {
                    if (Select_G | Select_Z) {
                      cat(sprintf("iteration %d: penalized log-likelihood = %.3f (observed-data = %.3f)\n",
                                  itr, new.loglik, raw.loglik))
                    }
                    else {
                      cat(sprintf("iteration %d: log-likelihood = %.3f\n",
                                  itr, new.loglik))
                    }
                  }
                  loglik_trace <- c(loglik_trace, raw.loglik)
                  # The I-step is a majorization step, not a true M-step, so
                  # under sporadic missingness the observed log-likelihood can
                  # dip by a small, bounded amount without indicating a real
                  # problem (measured worst case ~-0.16 in the test suite);
                  # without imputation the ascent should hold to numerical
                  # tolerance. Anything past that slack is flagged.
                  # Selection-driven M-steps (glasso for Rho_Z_Cov, thresholding
                  # for Rho_G/Rho_Z_Mu) are not a guaranteed ascent step either,
                  # and empirically dip by comparable amounts to the I-step
                  # (observed up to ~-0.41 in the test suite) -- both regimes
                  # get the wider slack. Even without either, floating-point
                  # accumulation (e.g. a downstream serial stage whose "G" is a
                  # continuous posterior probability rather than a clean
                  # indicator) can produce dips of a few thousandths, so the
                  # "clean" slack is a small multiple of tol rather than tol
                  # itself.
                  dip_slack <- if (isTRUE(na_pattern$impute_flag) || penalty_requested) 0.5 else max(tol, 0.01)
                  if (is.finite(res.loglik) && is.finite(new.loglik) &&
                      new.loglik < res.loglik - dip_slack) {
                    warning(sprintf(
                      "LUCID early model: log-likelihood decreased at iteration %d (from %.6f to %.6f); this indicates a numerical problem, not expected EM behaviour.",
                      itr, res.loglik, new.loglik), call. = FALSE)
                  }
                  if (check_convergence(res.loglik, new.loglik, abs_tol = tol)) {
                    convergence <- TRUE
                    em_converged <- TRUE
                    em_iter <- itr
                    if (verbose) {
                      cat("Success: LUCID early model converged!",
                        "\n\n")
                    }
                  }
                  res.loglik <- new.loglik
                }
            }
            if (!convergence && !restart_needed) {
                warning("LUCID early model reached max_itr = ", max_itr,
                        " without meeting the convergence tolerance (tol = ", tol,
                        "); returning the latest estimates. Inspect fit$em_control$converged.",
                        call. = FALSE)
                em_converged <- FALSE
                em_iter <- itr
                convergence <- TRUE
            }
        }
        res.likelihood <- Estep_early(beta = res.beta, mu = res.mu, 
            sigma = res.sigma, gamma = res.gamma, G = G, Z = Z, 
            Y = Y, family.list = family.list, itr = itr, CoY = CoY, 
            N = N, K = K, dimCoY = dimCoY, useY = useY, ind.na = na_pattern$indicator_na)
        res.r <- t(apply(res.likelihood, 1, lse_vec))
        if (!useY) {
            res.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY,
                K = K, CoYnames = CoYnames)
        }
        res.loglik <- observed_loglik(res.likelihood)
        pars <- if (isTRUE(relabel_by_outcome)) {
            switch_Y(beta = res.beta, mu = res.mu, sigma = res.sigma,
                     gamma = res.gamma, K = K)
        } else {
            # Unsupervised stage: order deterministically by the omics cluster
            # means so the labels (and every coefficient reported against them)
            # are reproducible across seeds.  See omics_order_index().
            relabel_early_parameters(beta = res.beta, mu = res.mu,
                                     sigma = res.sigma, gamma = res.gamma, K = K,
                                     index = omics_order_index(res.mu))
        }
        res.r <- res.r[, pars$index]
        colnames(pars$beta) <- c("intercept", Gnames)
        colnames(pars$mu) <- Znames
        if (Select_G) {
            beta_exposure <- pars$beta[, 1 + seq_len(dimG), drop = FALSE]
            coef_range <- apply(beta_exposure, 2, function(x) diff(range(x)))
            selectG <- as.logical(abs(coef_range) > 0.001)
        }
        else {
            selectG <- rep(TRUE, dimG)
        }
        if (Select_Z) {
            tt2 <- apply(pars$mu, 2, range)
            selectZ <- abs(tt2[2, ] - tt2[1, ]) > 0.001
        }
        else {
            selectZ <- rep(TRUE, dimZ)
        }
        if (isTRUE(verbose)) {
            if (penalty_requested) {
                cat(sprintf("Finished LUCID early model: observed-data log-likelihood = %.3f; selected G = %d/%d; selected Z = %d/%d.\n\n",
                    res.loglik, sum(selectG), length(selectG), sum(selectZ), length(selectZ)))
            }
            else {
                cat(sprintf("Finished LUCID early model: log-likelihood = %.3f.\n\n",
                    res.loglik))
            }
        }
        else {
            if (penalty_requested) {
                cat(sprintf("Finished LUCID early model. Selected G: %d/%d; Selected Z: %d/%d.\n", 
                    sum(selectG), length(selectG), sum(selectZ), length(selectZ)))
            }
            else {
                cat("Finished LUCID early model.\n")
            }
        }
        results <- list(res_Beta = pars$beta, res_Mu = pars$mu, 
            res_Sigma = pars$sigma, res_Gamma = pars$gamma, K = K, 
            var.names = list(Gnames = Gnames, Znames = Znames, 
                Ynames = Ynames), init_omic.data.model = model.best, 
            likelihood = res.loglik, inclusion.p = res.r, family = family, 
            select = list(selectG = selectG, selectZ = selectZ), 
            useY = useY, Z = Z, init_impute = init_impute, init_par = init_par, 
            Rho = list(Rho_G = Rho_G, Rho_Z_Mu = Rho_Z_Mu, Rho_Z_Cov = Rho_Z_Cov), 
            missing_summary = missing_summary,
            em_control = list(tol = tol, max_itr = max_itr, max_tot.itr = max_tot.itr,
                converged = em_converged, n_iter = em_iter,
                n_restart = n_restart, loglik_trace = loglik_trace,
                n_starts = 1L, n_starts_ok = 1L))
        class(results) <- c("early_lucid")
        return(results)
    }
    else {
        if (is.null(Z)) {
            stop("Input data 'Z' is missing")
        }
        if (!is.list(Z)) {
            stop("Input data 'Z' should be a list for LUCID in Parallel!")
        }
        else {
            for (i in 1:length(Z)) {
                if (!is.matrix(Z[[i]])) {
                  Z[[i]] <- as.matrix(Z[[i]])
                  if (!is.numeric(Z[[i]])) {
                    stop("Input data 'Z' should be numeric")
                  }
                }
            }
        }
        Znames <- vector("list", length(Z))
        for (i in 1:length(Z)) {
            if (is.null(colnames(Z[[i]]))) {
                Znames[[i]] <- paste0("Z_", i, "_", 1:ncol(Z[[i]]))
            }
            else {
                Znames[[i]] <- colnames(Z[[i]])
            }
        }
        N <- nrow(G)
        nOmics <- length(Z)
        nG <- ncol(G)
        nZ <- as.integer(sapply(Z, ncol))
        family <- to_parallel_family(family)
        # NULL means "let mclust choose", and is resolved per layer at
        # initialization below. A single name applies to every layer; a
        # per-layer vector is recycled to exactly one entry per layer
        # (`rep()` would have squared its length).
        modelNames <- if (is.null(init_omic.data.model)) {
            NULL
        } else {
            rep_len(init_omic.data.model, length(K))
        }
        dimCoG <- ifelse(is.null(CoG), 0, ncol(CoG))
        dimCoY <- ifelse(is.null(CoY), 0, ncol(CoY))
        if (!is.null(CoG)) {
            G <- cbind(G, CoG)
            Gnames <- c(Gnames, CoGnames)
        }
        na_pattern <- vector("list", nOmics)
        for (i in 1:nOmics) {
            na_pattern[[i]] <- check_na(Z[[i]])
            if (na_pattern[[i]]$impute_flag) {
                if (init_impute == "mix") {
                  if (!requireNamespace("mix", quietly = TRUE)) {
                    stop("init_impute = \"mix\" requires the 'mix' package, ",
                      "which is archived on CRAN; install it manually (e.g. ",
                      "from the CRAN Archive) or use init_impute = \"lod\" ",
                      "instead.", call. = FALSE)
                  }
                  if (verbose) {
                    cat("Intializing imputation of missing values in 'Z' via the mix package \n\n")
                  }
                  invisible(capture.output(Z[[i]] <- mclust::imputeData(Z[[i]],
                    seed = seed)))
                  Z[[i]][na_pattern[[i]]$indicator_na == 3, ] <- NA
                }
                if (init_impute == "lod") {
                  if (verbose) {
                    cat("Intializing imputation of missing values in 'Z' via LOD / sqrt(2) \n\n")
                  }
                  Z[[i]] <- apply(Z[[i]], 2, fill_data_lod)
                  colnames(Z[[i]]) <- Znames[[i]]
                }
                # D7: see the early branch
                Z[[i]][na_pattern[[i]]$indicator_na == 3, ] <- NA
            }
        }
        missing_summary <- summarize_missing_stats(list(index = lapply(na_pattern, `[[`, "index"), 
            indicator_na = lapply(na_pattern, `[[`, "indicator_na")), lucid_model = "parallel")
        if (!isTRUE(verbose)) {
            cat(sprintf("Fitting LUCID parallel model (%d layers)...\n", nOmics))
        }
        tot.itr <- 0
        flag_converge <- FALSE
        em_converged <- FALSE
        em_iter <- 0L
        loglik_trace <- numeric(0)
        n_restart <- 0L
        while (!flag_converge && tot.itr <= max_tot.itr) {
            if (tot.itr > 0) {
                n_restart <- n_restart + 1L
            }
            if (tot.itr > 0) {
                seed <- seed + 10
            }
            set.seed(seed)
            if (init_par == "mclust") {
                Mu_Sigma <- initialize_Mu_Sigma(K = K, Z = Z, modelNames = modelNames,
                  na_pattern = na_pattern)
                Mu <- Mu_Sigma$Mu
                Sigma <- Mu_Sigma$Sigma
                init_z <- Mu_Sigma$z
                # Adopt whatever mclust selected per layer. Without this the
                # M-step below receives modelNames = NULL and mclust::mstep()
                # fails on `switch(EXPR = modelName, ...)`, so the documented
                # automatic selection was unusable for the parallel model.
                modelNames <- Mu_Sigma$modelNames
            }
            else {
                if (is.null(modelNames)) {
                  # Random initialization has no mclust fit to select from, so
                  # fall back to the same default the early path uses.
                  modelNames <- rep("EEV", nOmics)
                  if (verbose) {
                    cat("GMM model for LUCID is not specified, 'EEV' model is used by default \n")
                  }
                }
                Mu <- initialize_Mu(K = K, nZ = nZ)
                Sigma <- initialize_Sigma(K = K, nZ = nZ)
                init_z <- vector(mode = "list", length = nOmics)
                for (i in 1:nOmics) {
                  temp_z <- matrix(runif(N * K[i]), nrow = N, ncol = K[i])
                  init_z[[i]] <- t(apply(temp_z, 1, function(x) x/sum(x)))
                }
            }
            Beta <- vector(mode = "list", length = nOmics)
            for (i in 1:nOmics) {
                invisible(capture.output(temp_fit <- nnet::multinom(init_z[[i]] ~ 
                  G)))
                Beta[[i]] <- coef(temp_fit)
            }
            Gamma <- initialize_Delta(K = K, CoY = CoY, family = family, 
                z = init_z, Y = Y)
            loglik <- -Inf
            if (verbose) {
                cat("Fitting LUCID in Parallel model", paste0("(", 
                  "K = ", K, ", Rho_G = ", Rho_G, ", Rho_Z_Mu = ", 
                  Rho_Z_Mu, ", Rho_Z_Cov = ", Rho_Z_Cov, ")"), 
                  "\n")
            }
            itr <- 0
            loglik_update <- -Inf
            restart_needed <- FALSE
            while (!flag_converge & itr < max_itr) {
                itr <- itr + 1
                tot.itr <- tot.itr + 1
                Estep_array <- Estep(G = G, Z = Z, Y = Y, Beta = Beta, 
                  Mu = Mu, Sigma = Sigma, Delta = Gamma, family = family, 
                  useY = useY, na_pattern = na_pattern, CoY = CoY)
                Estep_r <- Estep_to_r(Estep_array = Estep_array, 
                  K = K, N = N)
                if (!all(is.finite(Estep_r))) {
                  if (verbose) {
                    cat("iteration", itr, ": EM algorithm collapsed: invalid estiamtes due to over/underflow, try LUCID with another seed \n")
                  }
                  restart_needed <- TRUE
                  break
                }
                else {
                  if (isTRUE(verbose)) {
                    cat("iteration", itr, ": E-step finished.\n")
                  }
                }
                res_Beta <- Mstep_GtoX(G = G, r = Estep_r, selectG = Select_G, 
                  penalty = Rho_G, K = K, N = N, dimG_exposure = nG)
                res_Mu_Sigma <- Mstep_XtoZ(Z = Z, r = Estep_r, 
                  K = K, modelNames = modelNames, N = N, na_pattern = na_pattern, 
                  selectZ = Select_Z, penalty.mu = Rho_Z_Mu, penalty.cov = Rho_Z_Cov,
                  mu = Mu)
                if (useY) {
                  res_Gamma <- Mstep_XtoY(Y = Y, CoY = CoY, r = Estep_r, 
                    K = K, N = N, family = family)
                }
                if (is.null(res_Mu_Sigma$Mu)) {
                  if (verbose) {
                    cat("variable selection failed, try LUCID with another seed \n")
                  }
                  restart_needed <- TRUE
                  break
                }
                  check.value <- all(is.finite(unlist(res_Beta$Beta)),
                  is.finite(unlist(res_Mu_Sigma$Mu)), if (useY) {
                    is.finite(c(res_Gamma$Gamma$coef, res_Gamma$Gamma$sd))
                  })
                if (!check.value) {
                  if (verbose) {
                    cat("iteration", itr, ": Invalid estimates, try LUCID with another seed \n")
                  }
                  restart_needed <- TRUE
                  break
                }
                else {
                  Beta <- res_Beta$Beta
                  Mu <- res_Mu_Sigma$Mu
                  Sigma <- res_Mu_Sigma$Sigma
                  if (useY) {
                    Gamma <- res_Gamma$Gamma
                  }
                }
                # Recompute responsibilities at the UPDATED parameters before
                # the I-step, for the same ascent reason as the early path.
                istep_array <- try(Estep(G = G, Z = Z, Y = Y, Beta = Beta,
                  Mu = Mu, Sigma = Sigma, Delta = Gamma, family = family,
                  useY = useY, na_pattern = na_pattern, CoY = CoY), silent = TRUE)
                istep_r <- if (inherits(istep_array, "try-error")) {
                  Estep_r
                } else {
                  tmp <- Estep_to_r(istep_array, K = K, N = N)
                  if (all(is.finite(tmp))) tmp else Estep_r
                }
                post.p <- vector(mode = "list", length = nOmics)
                for (i in 1:nOmics) {
                  post.p[[i]] = compute_res_r(r = istep_r, N = N, 
                    layer = i)
                }
                for (i in 1:nOmics) {
                  if (na_pattern[[i]]$impute_flag) {
                    Z[[i]] <- Istep_Z(Z = Z[[i]], p = post.p[[i]], 
                      mu = Mu[[i]], sigma = Sigma[[i]], index = na_pattern[[i]]$index, 
                      lucid_model = "parallel")
                  }
                }
                updated_array <- Estep(G = G, Z = Z, Y = Y, Beta = Beta,
                  Mu = Mu, Sigma = Sigma, Delta = Gamma, family = family,
                  useY = useY, na_pattern = na_pattern, CoY = CoY)
                Estep_r <- Estep_to_r(updated_array, K = K, N = N)
                loglik_update <- cal_loglik(Estep_array = updated_array)
                if (isTRUE(verbose)) {
                  cat(sprintf("iteration %d: log-likelihood = %.3f\n", itr, loglik_update))
                }
                loglik_trace <- c(loglik_trace, loglik_update)
                # Same rationale as the early path: the I-step is a
                # majorization step, so a bounded dip under sporadic
                # missingness is expected (measured worst case ~-0.04 in the
                # test suite) and only a decrease beyond that slack is flagged.
                # Same rationale as the early path (see comment there): a
                # selection-driven M-step is not a guaranteed ascent either,
                # so it gets the same wider slack as the I-step regime, and the
                # "clean" slack allows for ordinary floating-point noise.
                dip_slack <- if (any(vapply(na_pattern, function(np) isTRUE(np$impute_flag), logical(1))) ||
                                  penalty_requested) {
                  0.5
                } else {
                  max(tol, 0.01)
                }
                if (is.finite(loglik) && is.finite(loglik_update) &&
                    loglik_update < loglik - dip_slack) {
                  warning(sprintf(
                    "LUCID parallel model: log-likelihood decreased at iteration %d (from %.6f to %.6f); this indicates a numerical problem, not expected EM behaviour.",
                    itr, loglik, loglik_update), call. = FALSE)
                }
                if (check_convergence(loglik, loglik_update, abs_tol = tol)) {
                  flag_converge <- TRUE
                  em_converged <- TRUE
                  em_iter <- itr
                  if (verbose) {
                    cat("Success: LUCID parallel model converged!",
                      "\n\n")
                  }
                }
                else {
                  loglik <- loglik_update
                }
            }
            if (!flag_converge && !restart_needed) {
                warning("LUCID parallel model reached max_itr = ", max_itr,
                        " without meeting the convergence tolerance (tol = ", tol,
                        "); returning the latest estimates. Inspect fit$em_control$converged.",
                        call. = FALSE)
                em_converged <- FALSE
                em_iter <- itr
                flag_converge <- TRUE
            }
        }
        selectG_layer <- res_Beta$selectG
        if (is.list(selectG_layer) && length(selectG_layer) > 0) {
            selectG <- Reduce("|", selectG_layer)
        }
        else {
            selectG <- selectG_layer
        }
        selectZ <- res_Mu_Sigma$selectZ
        if (!useY) {
            Estep_array <- Estep(G = G, Z = Z, Y = Y, Beta = Beta,
                Mu = Mu, Sigma = Sigma, Delta = Gamma, family = family,
                useY = FALSE, na_pattern = na_pattern, CoY = CoY)
            Estep_r <- Estep_to_r(Estep_array = Estep_array, K = K, N = N)
            res_Gamma <- Mstep_XtoY(Y = Y, CoY = CoY, r = Estep_r, 
                K = K, N = N, family = family)
            Gamma <- res_Gamma$Gamma
            loglik_update <- cal_loglik(Estep_array = Estep_array)
        }
        relabeled <- if (isTRUE(relabel_by_outcome)) {
            relabel_parallel_parameters(Beta = Beta, Mu = Mu, Sigma = Sigma,
                Delta = Gamma, r = Estep_r, K = K, selectZ = selectZ)
        } else {
            # Unsupervised stage: deterministic ordering by omics cluster means
            # (see the early branch above).
            relabel_parallel_parameters(Beta = Beta, Mu = Mu, Sigma = Sigma,
                Delta = Gamma, r = Estep_r, K = K, selectZ = selectZ,
                permutations = lapply(Mu, omics_order_index))
        }
        Beta <- relabeled$Beta
        res_Beta$Beta <- Beta
        Mu <- relabeled$Mu
        Sigma <- relabeled$Sigma
        Gamma <- relabeled$Delta
        Estep_r <- relabeled$r
        selectZ <- relabeled$selectZ
        res_Gamma <- list(fit = Gamma$fit, Gamma = Gamma)
        post.p <- lapply(seq_along(K), function(i) compute_res_r(Estep_r, N, i))
        final_loglik <- ifelse(is.finite(loglik_update), loglik_update, 
            loglik)
        if (isTRUE(verbose)) {
            if (penalty_requested) {
                z_sel <- vapply(selectZ, function(x) {
                  if (is.null(dim(x))) sum(x) else sum(colSums(x) > 0)
                }, numeric(1))
                z_tot <- vapply(selectZ, function(x) {
                  if (is.null(dim(x))) length(x) else ncol(x)
                }, numeric(1))
                cat(sprintf("Finished LUCID parallel model: observed-data log-likelihood = %.3f; selected G = %d/%d; selected Z by layer = %s.\n\n",
                    final_loglik, sum(selectG), length(selectG), paste0(z_sel, "/", z_tot, collapse = ", ")))
            }
            else {
                cat(sprintf("Finished LUCID parallel model: log-likelihood = %.3f.\n\n", 
                    final_loglik))
            }
        }
        else {
            if (penalty_requested) {
                z_sel <- vapply(selectZ, function(x) {
                  if (is.null(dim(x))) sum(x) else sum(colSums(x) > 0)
                }, numeric(1))
                z_tot <- vapply(selectZ, function(x) {
                  if (is.null(dim(x))) length(x) else ncol(x)
                }, numeric(1))
                cat(sprintf("Finished LUCID parallel model. Selected G: %d/%d; Selected Z by layer: %s.\n", 
                    sum(selectG), length(selectG), paste0(z_sel, "/", z_tot, collapse = ", ")))
            }
            else {
                cat("Finished LUCID parallel model.\n")
            }
        }
        results <- list(res_Beta = res_Beta, res_Mu = Mu, res_Sigma = Sigma, 
            res_Gamma = res_Gamma, K = K, N = N, var.names = list(Gnames = Gnames, 
                Znames = Znames, Ynames = Ynames), init_omic.data.model = modelNames, 
            likelihood = final_loglik, inclusion.p = post.p, 
            family = family, select = list(selectG = selectG, 
                selectG_layer = selectG_layer, selectZ = selectZ), useY = useY, Z = Z, z = Estep_r, 
            init_impute = init_impute, init_par = init_par, Rho = list(Rho_G = Rho_G, 
                Rho_Z_Mu = Rho_Z_Mu, Rho_Z_Cov = Rho_Z_Cov), missing_summary = missing_summary,
            em_control = list(tol = tol, max_itr = max_itr, max_tot.itr = max_tot.itr,
                converged = em_converged, n_iter = em_iter,
                n_restart = n_restart, loglik_trace = loglik_trace,
                n_starts = 1L, n_starts_ok = 1L))
        class(results) <- c("lucid_parallel")
        return(results)
    }
}

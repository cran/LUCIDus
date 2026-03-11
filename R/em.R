############workhorse function for estimate_lucid############
############EM algorithm for "early", "parallel"############
est_lucid <- 
function (lucid_model = c("early", "parallel"), G, Z, Y, CoG = NULL, 
    CoY = NULL, K, init_omic.data.model = "EEV", useY = TRUE, 
    tol = 0.001, max_itr = 1000, max_tot.itr = 10000, Rho_G = 0, 
    Rho_Z_Mu = 0, Rho_Z_Cov = 0, family = c("normal", "binary"), 
    seed = 123, init_impute = c("mix", "lod"), init_par = c("mclust", 
        "random"), verbose = FALSE) 
{
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
        }
        tot.itr <- 0
        convergence <- FALSE
        while (!convergence && tot.itr <= max_tot.itr) {
            if (tot.itr > 0) {
                seed <- seed + 10
            }
            set.seed(seed)
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
                  check.gamma <- is.finite(unlist(new.gamma))
                }
                if (na_pattern$impute_flag) {
                  Z <- Istep_Z(Z = Z, p = res.r, mu = res.mu, 
                    sigma = res.sigma, index = na_pattern$index, 
                    lucid_model = "early")
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
                  new.loglik <- sum(rowSums(res.r * new.likelihood))
                  if (Select_G) {
                    beta_exposure <- res.beta[, 1 + seq_len(dimG), drop = FALSE]
                    new.loglik <- new.loglik - Rho_G * sum(abs(beta_exposure))
                  }
                  if (Select_Z) {
                    new.loglik <- new.loglik - Rho_Z_Mu * sum(abs(res.mu)) - 
                      Rho_Z_Cov * sum(abs(res.sigma))
                  }
                  if (isTRUE(verbose)) {
                    if (Select_G | Select_Z) {
                      cat("iteration", itr, ": M-step finished, ", 
                        "penalized loglike = ", sprintf("%.3f", 
                          new.loglik), "\n")
                    }
                    else {
                      cat("iteration", itr, ": M-step finished, ", 
                        "loglike = ", sprintf("%.3f", new.loglik), 
                        "\n")
                    }
                  }
                  if (abs(res.loglik - new.loglik) < tol) {
                    convergence <- TRUE
                    if (verbose) {
                      cat("Success: Early Integration LUCID converges!", 
                        "\n\n")
                    }
                  }
                  res.loglik <- new.loglik
                }
            }
            if (!convergence && !restart_needed) {
                if (verbose) {
                  cat("Maximum EM iterations reached without convergence; returning latest estimates.\n")
                }
                convergence <- TRUE
            }
        }
        if (!useY) {
            res.gamma <- Mstep_Y(Y = Y, r = res.r, CoY = CoY, 
                K = K, CoYnames = CoYnames)
        }
        res.likelihood <- Estep_early(beta = res.beta, mu = res.mu, 
            sigma = res.sigma, gamma = res.gamma, G = G, Z = Z, 
            Y = Y, family.list = family.list, itr = itr, CoY = CoY, 
            N = N, K = K, dimCoY = dimCoY, useY = useY, ind.na = na_pattern$indicator_na)
        res.r <- t(apply(res.likelihood, 1, lse_vec))
        res.loglik <- sum(rowSums(res.r * res.likelihood))
        if (Select_G) {
            beta_exposure <- res.beta[, 1 + seq_len(dimG), drop = FALSE]
            res.loglik <- res.loglik - Rho_G * sum(abs(beta_exposure))
        }
        if (Select_Z) {
            res.loglik <- res.loglik - Rho_Z_Mu * sum(abs(res.mu)) - 
                Rho_Z_Cov * sum(abs(res.sigma))
        }
        pars <- switch_Y(beta = res.beta, mu = res.mu, sigma = res.sigma, 
            gamma = res.gamma, K = K)
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
                cat(sprintf("Finished LUCID early model: penalized log-likelihood = %.3f; selected G = %d/%d; selected Z = %d/%d.\n\n", 
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
            em_control = list(tol = tol, max_itr = max_itr, max_tot.itr = max_tot.itr))
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
        modelNames = rep(init_omic.data.model, length(K))
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
            }
        }
        missing_summary <- summarize_missing_stats(list(index = lapply(na_pattern, `[[`, "index"), 
            indicator_na = lapply(na_pattern, `[[`, "indicator_na")), lucid_model = "parallel")
        if (!isTRUE(verbose)) {
            cat(sprintf("Fitting LUCID parallel model (%d layers)...\n", nOmics))
        }
        tot.itr <- 0
        flag_converge <- FALSE
        while (!flag_converge && tot.itr <= max_tot.itr) {
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
            }
            else {
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
                  useY = useY, na_pattern = na_pattern)
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
                  selectZ = Select_Z, penalty.mu = Rho_Z_Mu, penalty.cov = Rho_Z_Cov)
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
                    is.finite(unlist(res_Gamma$Gamma))
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
                post.p <- vector(mode = "list", length = nOmics)
                for (i in 1:nOmics) {
                  post.p[[i]] = compute_res_r(r = Estep_r, N = N, 
                    layer = i)
                }
                for (i in 1:nOmics) {
                  if (na_pattern[[i]]$impute_flag) {
                    Z[[i]] <- Istep_Z(Z = Z[[i]], p = post.p[[i]], 
                      mu = Mu[[i]], sigma = Sigma[[i]], index = na_pattern[[i]]$index, 
                      lucid_model = "parallel")
                  }
                }
                loglik_update <- cal_loglik(Estep_array = Estep_array, 
                  Estep_r = Estep_r)
                if (isTRUE(verbose)) {
                  cat(sprintf("iteration %d: log-likelihood = %.3f\n", itr, loglik_update))
                }
                if (abs(loglik - loglik_update) < tol) {
                  flag_converge <- TRUE
                  if (verbose) {
                    cat("Success: LUCID in parallel converges!", 
                      "\n\n")
                  }
                }
                else {
                  loglik <- loglik_update
                }
            }
            if (!flag_converge && !restart_needed) {
                if (verbose) {
                  cat("Maximum EM iterations reached without convergence; returning latest estimates.\n")
                }
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
            res_Gamma <- Mstep_XtoY(Y = Y, CoY = CoY, r = Estep_r, 
                K = K, N = N, family = family)
            Gamma <- res_Gamma$Gamma
        }
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
                cat(sprintf("Finished LUCID parallel model: penalized log-likelihood = %.3f; selected G = %d/%d; selected Z by layer = %s.\n\n", 
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
            em_control = list(tol = tol, max_itr = max_itr, max_tot.itr = max_tot.itr))
        class(results) <- c("lucid_parallel")
        return(results)
    }
}

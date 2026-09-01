#' @title Predict Cluster Assignment and Outcome From a Fitted LUCID Model
#' @description Predict cluster assignment and outcome using new data on G, Z, and optional Y.
#' If \code{g_computation = TRUE}, prediction uses only the G-to-X path from the
#' fitted model and returns counterfactual-style predictions under modified G.
#' This function can also be used to extract latent cluster assignments when using
#' the training data as input.
#' @param model A model fitted and returned by \code{\link{estimate_lucid}}
#' @param lucid_model Optional; "early", "parallel", or "serial". Auto-detected
#' from \code{class(model)} when omitted (the normal case), so this rarely
#' needs to be set explicitly -- it exists for backward compatibility with
#' scripts written before auto-detection. A serial model must have at least
#' two stages to be predicted; a single-stage serial model is a fully
#' equivalent early or parallel model and should be fitted as one.
#' @param G Exposures, a numeric vector, matrix, or data frame. Categorical variable
#' should be transformed into dummy variables. If a matrix or data frame, rows
#' represent observations and columns correspond to variables.
#' @param Z Omics data, and **required for every model type** unless
#' \code{g_computation = TRUE}. If "early", an N by M matrix. If "parallel", a
#' list, each element i is a matrix with N rows and P_i features. If "serial", a
#' list, each element i is a matrix with N rows and p_i features (or a list with
#' two or more matrices with N rows and a certain number of features).
#'
#' The requirement is not arbitrary: the E-step forms the posterior from the
#' omics likelihood, so with no \code{Z} there is nothing to condition on.
#' \code{g_computation = TRUE} is a different estimator, not a way around this
#' -- it drops the omics and outcome terms and uses the exposure path alone --
#' which is why it is the one mode that accepts \code{Z = NULL}.
#' @param Y Outcome, a numeric vector. Categorical variable is not allowed. Binary
#' outcome should be coded as 0 and 1.
#' @param CoG Optional, covariates to be adjusted for estimating the latent cluster.
#' A numeric vector, matrix or data frame. Categorical variable should be transformed
#' into dummy variables.
#' @param CoY Optional, covariates to be adjusted for estimating the association
#' between latent cluster and the outcome. A numeric vector, matrix or data frame.
#' Categorical variable should be transformed into dummy variables.
#' @param response If \code{TRUE}, when predicting binary outcomes, class labels
#' (0/1) are returned using a 0.5 threshold. If \code{FALSE}, predicted
#' probabilities are returned.
#' @param g_computation If \code{TRUE}, prediction uses only information on G,
#' making it the counterfactual mode: hold the fitted model fixed, vary
#' \code{G}, and read off what the model implies. It is the only mode in which
#' \code{Z} may be omitted, and it is also the only one that returns
#' \code{pred.z}. Supplied \code{Z} and \code{Y} are ignored (with a printed
#' notice) for "early", "parallel", and "serial", so results are unchanged by
#' passing them.
#' @param verbose A flag indicates whether detailed information
#' is printed in console. Default is FALSE. Applies consistently to all three
#' model types (early, parallel, serial).
#' @return A list containing:
#' \item{inclusion.p}{Posterior inclusion probabilities for latent clusters (a
#' matrix for "early"; a list by layer for "parallel" and "serial"). Columns are
#' ordered by cluster, matching the row order of the model's \code{res_Mu} and
#' \code{res_Beta}.}
#' \item{pred.x}{Predicted latent-cluster labels (a numeric vector for "early";
#' a list by layer for "parallel" and "serial"), obtained as the maximum a
#' posteriori column of \code{inclusion.p}. Labels run \code{1, ..., K},
#' agreeing with Eq 21 and with the cluster names used by \code{summary()} and
#' the \code{mu} and \code{beta} row names. Versions before 3.1.0 returned
#' \code{0, ..., K - 1} here; code that compensated by adding one must drop
#' that adjustment.}
#' \item{pred.y}{Predicted outcome values. For binary outcomes, this is class
#' labels when \code{response = TRUE} and probabilities when
#' \code{response = FALSE}.}
#' \item{pred.z}{Predicted omics means under g-computation mode
#' (\code{g_computation = TRUE}); \code{NULL} otherwise.}
#'
#' Supplying \code{Y} makes the cluster prediction supervised: the outcome
#' enters the posterior alongside G and Z, as it does during fitting. Omitting
#' it predicts clusters from G and Z alone, which is what is wanted when the
#' outcome is unavailable or must not inform the assignment.
#' 
#' @export
#'
#' @examples
#' # prepare data (a small subset keeps the example quick)
#' G <- sim_data$G[1:150, ]
#' Z <- sim_data$Z[1:150, ]
#' Y_normal <- sim_data$Y_normal[1:150, , drop = FALSE]
#'
#' # fit lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", K = 2,
#'                        family = "normal", max_itr = 20, max_tot.itr = 50)
#'
#' # prediction on training set (lucid_model is auto-detected from fit1's class)
#' pred1 <- predict_lucid(model = fit1, G = G, Z = Z, Y = Y_normal)
#' pred2 <- predict_lucid(model = fit1, G = G, Z = Z)
#'
#' # g-computation style prediction using only G
#' pred_g <- predict_lucid(model = fit1, G = G, Z = NULL, g_computation = TRUE)
#'


predict_lucid <- function(model,
                          lucid_model = NULL,
                          G,
                          Z = NULL,
                          Y = NULL,
                          CoG = NULL,
                          CoY = NULL,
                          response = TRUE,
                          g_computation = FALSE,
                          verbose = FALSE){

  # Resolve the model type ONCE, here. `model`'s own class already says
  # whether it's early/parallel/serial, so lucid_model is auto-detected from
  # it by default; a caller may still name it explicitly (e.g. for backward
  # compatibility with older scripts), in which case it is validated against
  # the usual three choices, matching the prior behavior exactly.
  lucid_model <- if (is.null(lucid_model)) {
    .detect_lucid_model(model)
  } else {
    match.arg(lucid_model, c("early", "parallel", "serial"))
  }

  if (g_computation == TRUE){
    if (verbose && (!is.null(Z) || !is.null(Y))){
      cat("G-computation only uses input for G, and the G-to-X association, input of Z and Y will not be used for prediction.\n")
    }
  }

  # One rule for every model type: the omics data is required. The ordinary
  # E-step forms the posterior from the omics likelihood, so there is nothing to
  # condition on without it; g-computation is a different estimator that uses the
  # exposure path alone, and is the only mode that may omit Z.
  #
  # This check sits before the model-specific dispatch so all three types give
  # the same message. The serial branch used to compare length(Z) against
  # length(K) first, and length(NULL) is 0, so a missing Z was reported as a
  # topology mismatch and the serial is.null() branch below was unreachable.
  if (!isTRUE(g_computation) && is.null(Z)) {
    stop("Input data 'Z' is required for prediction. ",
         "Omit it only with g_computation = TRUE, which predicts from the ",
         "exposures alone.", call. = FALSE)
  }
  
  if (lucid_model == "early" | lucid_model == "parallel"){
    # ========================== Early Integration ==========================
    # ========================== LUCID IN PARALLEL ==========================
    res_pred = pred_lucid(model = model, lucid_model = lucid_model, G = G, Z = Z, Y = Y,
                          CoG = CoG, CoY = CoY, response = response, g_computation = g_computation,
                          verbose = verbose)
    return(res_pred)
  }else if (lucid_model == "serial"){
    # ========================== LUCID IN Serial ==========================
    n <- nrow(G)
    K <- model$K

    # A one-stage serial model is degenerate -- it is an early or parallel model
    # with extra wrapping -- and the stage loop below cannot predict it: the
    # branches are first / middle / last, but a single stage is simultaneously
    # first and last, so it takes the first branch and never reaches the code
    # that assigns pred.y. Rather than let that surface as
    # "object 'pred.y' not found", decline it and say what to do instead.
    if (length(K) < 2L) {
      stop("A serial model needs at least two stages to predict. ",
           "For a single stage, fit with lucid_model = \"early\" or ",
           "\"parallel\" instead.", call. = FALSE)
    }

    ## check data format ==== special for Z under serial
    # For g-computation, Z is optional and ignored.
    if (!g_computation) {
      # Order matters: a missing Z is caught by the guard above, a non-list Z
      # must be reported as such before its length is compared against K, and
      # only then does a genuine topology mismatch make sense as a message.
      if(is.null(Z)) {
        stop("Input data 'Z' is missing")
      }
      if(!is.list(Z)) {
        stop("Input data 'Z' should be a list for LUCID in Serial!")
      }
      if(length(Z) != length(K)) {
        stop("Z and K should be two lists of the same length for LUCID in Serial!")
      }
      {
        for(i in 1:length(K)) {
          if(is.numeric(K[[i]])) {
            if(!is.matrix(Z[[i]])) {
              stop("For LUCID in Serial, input data 'Z' must match the K input. When the element of K is a integer, the corresponding element of Z must also be a matrix!")
            }}
          if(is.list(K[[i]])) {
            if(!is.list(Z[[i]])) {
              stop("For LUCID in Serial, input data 'Z' must match the K input. When the element of K is a list, the corresponding element of Z must also be a list of matrices!")
            }
          }
        }
      }
    }

    # initiate the empty lists to store the predictions for each sub model
    post.p.list <- vector(mode = "list", length = length (K))
    pred.x.list <- vector(mode = "list", length = length (K))
    pred.z.list <- vector(mode = "list", length = length (K))
    
    #loop through each K
    for (i in 1:length(K)){
      if(verbose){
        cat(sprintf("Predicting LUCID serial model (Stage %d/%d)\n", i, length(K)))
      }
      ##Scenario 1: the first serial sub model
      if (i == 1){
        if (is.numeric(K[[1]])){
          #if the first serial sub model is early integration (1 layer)
          temp_pred = pred_lucid(model = model$submodel[[1]], lucid_model = "early", G = G, Z = if (g_computation) NULL else Z[[1]], Y = NULL,
                                CoG = CoG, CoY = NULL, g_computation = g_computation, response = FALSE)

          post.p.list[[1]] = temp_pred$inclusion.p
          pred.x.list[[1]] = temp_pred$pred.x
          if (g_computation == TRUE){
            pred.z.list[[1]] = temp_pred$pred.z
          }
          post.p = temp_pred$inclusion.p[,-1]

        }else{
          #if the first serial sub model is lucid in parallel
          temp_pred = pred_lucid(model = model$submodel[[1]], lucid_model = "parallel", G = G, Z = if (g_computation) NULL else Z[[1]], Y = NULL,
                                 CoG = CoG, CoY = NULL, g_computation = g_computation, response = FALSE)

          post.p.list[[1]] = temp_pred$inclusion.p
          pred.x.list[[1]] = temp_pred$pred.x
          if (g_computation == TRUE){
            pred.z.list[[1]] = temp_pred$pred.z
          }
          temp.p = temp_pred$inclusion.p
          temp.p.list = vector(mode = "list", length = length(temp.p))
          for (j in 1:length(temp.p)){
            temp.p.list[[j]] = temp.p[[j]][,-1]
          }
          post.p = matrix(unlist(temp.p.list), nrow = nrow(G), byrow = FALSE)

        }
      }else if (i < length(K)){
        ##Scenario 2: the middle serial sub models
        if (is.numeric(K[[i]])){
          #if the middle serial sub model is early integration (1 layer)
          temp_pred = pred_lucid(model = model$submodel[[i]], lucid_model = "early", G = post.p, Z = if (g_computation) NULL else Z[[i]], Y = NULL,
                                 CoG = NULL, CoY = NULL, g_computation = g_computation, response = FALSE)

          post.p.list[[i]] = temp_pred$inclusion.p
          pred.x.list[[i]] = temp_pred$pred.x
          if (g_computation == TRUE){
            pred.z.list[[i]] = temp_pred$pred.z
          }
          post.p = temp_pred$inclusion.p[,-1]

        }else{
          #if the first serial sub model is lucid in parallel
          temp_pred = pred_lucid(model = model$submodel[[i]], lucid_model = "parallel", G = post.p, Z = if (g_computation) NULL else Z[[i]], Y = NULL,
                                 CoG = NULL, CoY = NULL, g_computation = g_computation, response = FALSE)

          post.p.list[[i]] = temp_pred$inclusion.p
          pred.x.list[[i]] = temp_pred$pred.x
          if (g_computation == TRUE){
            pred.z.list[[i]] = temp_pred$pred.z
          }
          
          temp.p = temp_pred$inclusion.p
          temp.p.list = vector(mode = "list", length = length(temp.p))
          for (j in 1:length(temp.p)){
            temp.p.list[[j]] = temp.p[[j]][,-1]
          }
          post.p = matrix(unlist(temp.p.list), nrow = nrow(G), byrow = FALSE)

        }
      }else if (i == length(K)){
        ##Scenario 3: the last sub model
        if (is.numeric(K[[i]])){
          #if the last serial sub model is early integration (1 layer)
          temp_pred = pred_lucid(model = model$submodel[[i]], lucid_model = "early", G = post.p, Z = if (g_computation) NULL else Z[[i]], Y = if (g_computation) NULL else Y,
                                 CoG = NULL, CoY = CoY, g_computation = g_computation, response = response)

          post.p.list[[i]] = temp_pred$inclusion.p
          pred.x.list[[i]] = temp_pred$pred.x
          if (g_computation == TRUE){
            pred.z.list[[i]] = temp_pred$pred.z
          }
          
          pred.y = temp_pred$pred.y

        }else{
          #if the last serial sub model is parallel (multiple layers)
          temp_pred = pred_lucid(model = model$submodel[[i]], lucid_model = "parallel", G = post.p, Z = if (g_computation) NULL else Z[[i]], Y = if (g_computation) NULL else Y,
                                 CoG = NULL, CoY = CoY, g_computation = g_computation, response = response)

          post.p.list[[i]] = temp_pred$inclusion.p
          pred.x.list[[i]] = temp_pred$pred.x
          if (g_computation == TRUE){
            pred.z.list[[i]] = temp_pred$pred.z
          }
          
          pred.y = temp_pred$pred.y

          }
        }
    }
    
    if (g_computation == FALSE){
      results <- list(inclusion.p = post.p.list,
                      pred.x = pred.x.list,
                      pred.y = pred.y)
    }else{
      results <- list(inclusion.p = post.p.list,
                      pred.x = pred.x.list,
                      pred.z = pred.z.list,
                      pred.y = pred.y)
    }
    
    return(results)
    }
  }


# =============================================================================
# pred_lucid(): the prediction workhorse for "early" and "parallel" (merged in
# from the former pred_lucid.R). predict_lucid() above is its only caller, for
# both top-level early/parallel prediction and each early/parallel leg of a
# serial model's stage-by-stage prediction.
# =============================================================================

#' Prediction workhorse for the "early" and "parallel" LUCID models
#'
#' Computes posterior cluster assignments and predicted outcomes for a
#' fitted early- or parallel-integration model. \code{predict_lucid()} is
#' this function's only caller, both for top-level early/parallel prediction
#' and for each early/parallel leg of a serial model's stage-by-stage
#' prediction.
#'
#' @param model A fitted \code{early_lucid} or \code{lucid_parallel} object.
#' @param lucid_model "early" or "parallel".
#' @param G Exposure data.
#' @param Z Omics data; required unless \code{g_computation = TRUE}.
#' @param Y Outcome data; optional (omitting it gives an unsupervised
#'   posterior).
#' @param CoG,CoY Optional covariates for the exposure and outcome models.
#' @param response Whether to return the outcome prediction on the response
#'   scale (vs. the linear predictor).
#' @param g_computation If \code{TRUE}, predicts from the exposure path alone
#'   (ignoring \code{Z}/\code{Y} if supplied) rather than the ordinary E-step.
#' @param verbose Whether to print progress messages.
#' @return A list with \code{inclusion.p} (posterior cluster probabilities),
#'   \code{pred.x} (predicted cluster), \code{pred.y} (predicted outcome),
#'   and, under g-computation, \code{pred.z} (predicted omics values).
#' @noRd
pred_lucid <- function(model,
                       lucid_model = c("early", "parallel"),
                       G,
                       Z = NULL,
                       Y = NULL,
                       CoG = NULL,
                       CoY = NULL,
                       response = TRUE,
                       g_computation = FALSE,
                       verbose = FALSE){
  model_family_std <- normalize_family_label(model$family)
  model_family_parallel <- to_parallel_family(model$family)
  ## 1.1 check data format ====
  if(is.null(G)) {
    stop("Input data 'G' is missing")
  } else {
    if(!is.matrix(G)) {
      G <- as.matrix(G)
      if(!is.numeric(G)) {
        stop("Input data 'G' should be numeric; categorical variables should be transformed into dummies")
      }
    }
  }

  CoGnames <- NULL
  if(!is.null(CoG)) {
    if(!is.matrix(CoG)) {
      CoG <- as.matrix(CoG)
      if(!is.numeric(CoG)) {
        stop("Input data 'CoG' should be numeric; categroical variables should be transformed into dummies")
      }
    }
    if(is.null(colnames(CoG))) {
      CoGnames <- paste0("CoG", 1:ncol(CoG))
    } else {
      CoGnames <- colnames(CoG)
    }
    colnames(CoG) <- CoGnames
  }

  CoYnames <- NULL
  if(!is.null(CoY)) {
    if(!is.matrix(CoY)) {
      CoY <- as.matrix(CoY)
      if(!is.numeric(CoY)) {
        stop("Input data 'CoY' should be numeric; categorical variables should be transformed into dummies")
      }
    }
    if(is.null(colnames(CoY))) {
      CoYnames <- paste0("CoY", 1:ncol(CoY))
    } else {
      CoYnames <- colnames(CoY)
    }
    colnames(CoY) <- CoYnames
  }

  if(!is.null(Y)) {
    if(!is.matrix(Y)) {
     Y <- as.matrix(Y)
      if(!is.numeric(Y)) {
        stop("Input data 'Y' should be numeric; binary outcome should be transformed them into dummies")
      }
      if(ncol(Y) > 1) {
        stop("Only continuous 'Y' or binary 'Y' is accepted")
      }
    }
    if(is_binary_family(model$family)) {
      if(!(all(Y %in% c(0, 1)))) {
        stop("Binary outcome should be coded as 0 and 1")
      }
    }
  }

  if (match.arg(lucid_model) == "early"){

  if(!inherits(model, "early_lucid")) {
    stop("model should be an object of early_lucid fitted by est_lucid")
  }

  if(verbose) {
    cat("Predicting LUCID early model...\n")
  }

  indicator_na <- NULL
  if (g_computation == FALSE) {
    if(is.null(Z)) {
      stop("Input data 'Z' is missing")
    } else {
      if(!is.matrix(Z)) {
        Z <- as.matrix(Z)
        if(!is.numeric(Z)) {
          stop("Input data 'Z' should be numeric")
        }
      }
    }
    if (ncol(Z) != ncol(model$res_Mu)) {
      stop("Input data 'Z' should have ", ncol(model$res_Mu),
           " columns (the number of omics features the model was fitted on), ",
           "but has ", ncol(Z), ".")
    }
    indicator_na <- check_na(Z)$indicator_na
  }


    n <- nrow(G)
    K <- model$K

    # model parameters
    beta <- model$res_Beta
    mu <- model$res_Mu
    Sigma <- model$res_Sigma
    Sigma.array <- array(as.numeric(unlist(Sigma)), dim = c(rep(ncol(mu), 2), K))
    gamma <- model$res_Gamma

    G <- cbind(G, CoG)
    dimCoY <- 0
    if(!is.null(CoY)){
      dimCoY <- ncol(CoY)
    }
    family.list <- switch(model_family_std,
                          normal = normal(K = K, dimCoY),
                          binary = binary(K = K, dimCoY))
    useY_flag <- ifelse(is.null(Y), FALSE, TRUE)
    # 1 - predict latent cluster
    if (g_computation == FALSE){
      res <- Estep_early(beta = beta,
                  mu = mu,
                  sigma = Sigma.array,
                  gamma = gamma,
                  G = G,
                  Z = Z,
                  Y = Y,
                  CoY = CoY,
                  family.list = family.list,
                  K = K,
                  N = n,
                  itr = 2,
                  dimCoY = dimCoY,
                  useY = useY_flag,
                  ind.na = indicator_na)
    }else{
      res <- Estep_early_g(beta = beta,
                         mu = mu,
                         sigma = Sigma.array,
                         gamma = gamma,
                         G = G,
                         Z = Z,
                         Y = Y,
                         CoY = CoY,
                         family.list = family.list,
                         K = K,
                         N = n,
                         itr = 2,
                         dimCoY = dimCoY,
                         useY = useY_flag,
                         ind.na = indicator_na)
    }

    # normalize the log-likelihood to probability
    res.r <- t(apply(res, 1, lse_vec))
    # predicted latent cluster
    # D6: cluster labels are 1..K, matching Eq 21, summary(), inclusion.p
    # columns and the mu/beta row names.  Earlier releases returned 0..K-1.
    pred.x <- sapply(1:n, function(x) return(nnet::which.is.max(res.r[x, ])))
    pred.x <- as.numeric(pred.x)
    cluster_eta <- matrix(rep(early_gamma_levels(gamma, K), each = n), nrow = n)
    covariate <- early_gamma_covariate(gamma, K)
    if (!is.null(CoY) && length(covariate)) {
      cluster_eta <- cluster_eta + as.numeric(as.matrix(CoY) %*% covariate)
    }
    cluster_response <- if (is_binary_family(model$family)) stats::plogis(cluster_eta) else cluster_eta
    pred.y <- rowSums(res.r * cluster_response)
    if (is_binary_family(model$family) && isTRUE(response)) pred.y <- as.numeric(pred.y > 0.5)
    
    if (g_computation == FALSE){
      results <- list(inclusion.p = res.r,
                      pred.x = pred.x,
                      pred.y = as.vector(pred.y))
    }else{
      #compute the weighted (by pred.x and estimated mu from model) predicted mu for each obs of each feature
      z_names <- colnames(mu)
      if (is.null(z_names)) {
        z_names <- paste0("Z", seq_len(ncol(mu)))
      }
      pred.z = matrix(NA, nrow = n, ncol = ncol(mu))
      colnames(pred.z) = z_names
      for (i in 1:n){
        p_i = res.r[i,]
        mu_i = colSums(p_i * mu)
        pred.z[i,] = mu_i
      }
      colnames(pred.z) = z_names
      
      results <- list(inclusion.p = res.r,
                      pred.x = pred.x,
                      pred.z = pred.z,
                      pred.y = as.vector(pred.y))
    }

    if(verbose) {
      cat("Finished predicting.\n")
    }
    return(results)

  }else if (match.arg(lucid_model) == "parallel"){
    ###LUCID in parallel#####
    if(!inherits(model, "lucid_parallel")) {
      stop("model should be an object of lucid_parallel fitted by est_lucid")
    }

    if(verbose) {
      cat("Predicting LUCID parallel model...\n")
    }

    if(is.null(colnames(G))){
      Gnames <- paste0("G", 1:ncol(G))
    } else {
      Gnames <- colnames(G)
    }
    colnames(G) <- Gnames


    N <- nrow(G)
    K <- model$K
    nOmics <- length(K)

    # combine G and CoG to adjust for CoG
    if(!is.null(CoG)) {
      G <- cbind(G, CoG)
      Gnames <- c(Gnames, CoGnames)
    }
    if (g_computation == FALSE) {
      if(is.null(Z)) {
        stop("Input data 'Z' is missing")
      }
      if(!is.list(Z)) {
        stop("Input data 'Z' should be a list for LUCID in Parallel!")
      }
      if (length(Z) != nOmics) {
        stop("Input data 'Z' should have the same number of layers as model$K for LUCID in Parallel!")
      }
      for(i in 1:length(Z)) {
        if(!is.matrix(Z[[i]])) {
          Z[[i]] <- as.matrix(Z[[i]])
          if(!is.numeric(Z[[i]])) {
            stop("Input data 'Z' should be numeric")
          }
        }
      }
    }

    #Get na_pattern, but for prediction, Z should be non-missing
    na_pattern <- vector("list", nOmics)
    if (g_computation == FALSE) {
      for(i in 1:nOmics) {
        na_pattern[[i]] <- check_na(Z[[i]])
      }
    }

    Beta = model$res_Beta$Beta
    Mu = model$res_Mu
    Sigma = model$res_Sigma
    Gamma <- model$res_Gamma$Gamma
    useY_flag <- ifelse(is.null(Y), FALSE, TRUE)
    
    if(g_computation == FALSE){
      Estep_array <- Estep(G = G, Z = Z, Y = Y,
                           Beta = Beta, Mu = Mu, Sigma = Sigma, Delta = Gamma,
                           family = model$family, useY = useY_flag, na_pattern = na_pattern,
                           CoY = CoY)

    }else{
      Estep_array <- Estep_g(G = G, Z = Z, Y = Y,
                           Beta = Beta, Mu = Mu, Sigma = Sigma, Delta = Gamma,
                           family = model$family, useY = useY_flag, na_pattern = na_pattern)
    }
    
    Estep_r <- Estep_to_r(Estep_array = Estep_array,
                          K = K,
                          N = N)
    
    post.p <- vector(mode = "list", length = nOmics)
    for(i in 1:nOmics) {
      post.p[[i]] = compute_res_r(r = Estep_r, N = N,layer = i)
    }

    # initialize container for predicted value
    pred_X <- vector(mode = "list", length = nOmics)
    # prediction of X and Y based on fitted data


      # 1 - prediction for X
      for (i in 1:nOmics) {
        # D6: labels are 1..K (see the early branch above)
        pred_X[[i]] <- sapply(1:N, function(x) return(nnet::which.is.max(post.p[[i]][x, ])))
      }

    # 2 - prediction for Y
    r = Estep_r

    state_eta <- parallel_state_eta(Gamma, CoY = CoY, N = N)
    state_response <- if(model_family_parallel == "binomial") stats::plogis(state_eta) else state_eta
    pred_Y <- colSums(matrix(r, nrow = prod(K), ncol = N) *
                      matrix(state_response, nrow = prod(K), ncol = N))
    if(model_family_parallel == "binomial" && isTRUE(response)) pred_Y <- as.numeric(pred_Y > 0.5)

    if (g_computation == FALSE){
      results <- list(inclusion.p = post.p,
                      pred.x = pred_X,
                      pred.y = pred_Y)
    }else{
      # Compute weighted predicted omics means using posterior cluster probabilities.
      # Handle both Mu storage conventions: K x Z and Z x K.
      pred.z <- vector(mode = "list", length = nOmics)
      for (i in 1:nOmics){
        mu_raw <- as.matrix(Mu[[i]])
        k_i <- ncol(post.p[[i]])
        if (nrow(mu_raw) == k_i) {
          mu_k_by_z <- mu_raw
          z_names <- colnames(mu_raw)
        } else if (ncol(mu_raw) == k_i) {
          mu_k_by_z <- t(mu_raw)
          z_names <- rownames(mu_raw)
        } else {
          stop("Incompatible dimensions between model$res_Mu and posterior inclusion probabilities in predict_lucid()")
        }

        z <- ncol(mu_k_by_z)
        if (is.null(z_names) || length(z_names) != z) {
          z_names <- paste0("Z_", i, "_", seq_len(z))
        }

        pred_layer_z <- matrix(NA_real_, nrow = N, ncol = z)
        colnames(pred_layer_z) <- z_names
        for (j in 1:N){
          p_j <- post.p[[i]][j, ]
          pred_layer_z[j, ] <- as.numeric(p_j %*% mu_k_by_z)
        }
        pred.z[[i]] <- pred_layer_z
      }
      results <- list(inclusion.p = post.p,
                      pred.x = pred_X,
                      pred.z = pred.z,
                      pred.y = pred_Y)
    }
      if(verbose) {
        cat("Finished predicting.\n")
      }
      return(results)
    }
  }

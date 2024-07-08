############workhorse function for predict_lucid############
############prediction for "early", "parallel"############

pred_lucid <- function(model,
                       lucid_model = c("early", "parallel"),
                       G,
                       Z,
                       Y = NULL,
                       CoG = NULL,
                       CoY = NULL,
                       response = TRUE,
                       g_computation = FALSE,
                       verbose = FALSE){
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
    if(model$family == "binary" | model$family == "binomial") {
      if(!(all(Y %in% c(0, 1)))) {
        stop("Binary outcome should be coded as 0 and 1")
      }
    }
  }

  if (match.arg(lucid_model) == "early"){

  if(!inherits(model, "early_lucid")) {
    stop("model should be an object of early_lucid fitted by est_lucid")
  }

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


    n <- nrow(G)
    K <- model$K

    # model parameters
    beta <- model$res_Beta
    mu <- model$res_Mu
    Sigma <- model$res_Sigma
    Sigma.array <- array(as.numeric(unlist(Sigma)), dim = c(rep(ncol(Z), 2), K))
    gamma <- model$res_Gamma

    G <- cbind(G, CoG)
    na_pattern <- check_na(Z)
    dimCoY <- 0
    if(!is.null(CoY)){
      dimCoY <- ncol(CoY)
    }
    family.list <- switch(model$family,
                          normal = normal(K = K, dimCoY),
                          binary = binary(K = K, dimCoY))
    useY_flag <- ifelse(is.null(Y), FALSE, TRUE)
    # browser()
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
                  ind.na = na_pattern$indicator_na)
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
                         ind.na = na_pattern$indicator_na)
    }

    # normalize the log-likelihood to probability
    res.r <- t(apply(res, 1, lse_vec))
    # predicted latent cluster
    pred.x <- sapply(1:n, function(x) return(nnet::which.is.max(res.r[x, ])))
    pred.x <- pred.x - 1
    #here, since we modify the gamma, where it's just intecept + betas instead of mu_Y in the model output
    #we need transfer gamma to mu_Y first 
    model_mu_Y <- gamma$beta
    model_mu_Y[2:K] <- model_mu_Y[1] + model_mu_Y[2:K]
    mu_Y <- cbind(res.r, CoY) %*% model_mu_Y
    #IP1 * beta0 + IP2 * (beta0 + beta1) + CoY*Gamma(CoY) is actually a weighed version of (Beta0 + beta1X)
    
    if(model$family == "normal"){
      #nomal Y is good
      pred.y <- mu_Y
    }
    if(model$family == "binary"){
      #need to discuss here with Dave, for normal, we derive mu_Y for each cluster and weighted by IP, 
      #for binary, do we get weight intercepts and betas and then map back to probability of the outcome??
      #Before, Yinqi use weighted_mu_logit to map back, 
      #For g_comp, how about keep logit of Y to estimate the causal effect of interest?
      
      #at first,I thought Yinqi is wrong
      #But now I think it makes sense, for exmaple if IP2 > IP1
      #then beta 0 is shrunk to 0
      #then IP1 * beta0 + IP2 * (beta0 + beta1) is actually a weighed version of (Beta0 + beta1X)
      
      #now I think it's good
      pred.y <- exp(mu_Y) / (1 + exp(mu_Y))
      if(response == TRUE){
        pred.y <- as.numeric(pred.y > 0.5)
      }
    }
    
    if (g_computation == FALSE){
      results <- list(inclusion.p = res.r,
                      pred.x = pred.x,
                      pred.y = as.vector(pred.y))
    }else{
      #compute the weighted (by pred.x and estimated mu from model) predicted mu for each obs of each feature
      z = ncol(mu)
      pred.z = matrix(NA, nrow = n, ncol = z)
      colnames(pred.z) = names(mu)
      for (i in 1:n){
        p_i = res.r[i,]
        mu_i = colSums(p_i * mu)
        pred.z[i,] = mu_i
      }
      colnames(pred.z) = names(mu)
      
      results <- list(inclusion.p = res.r,
                      pred.x = pred.x,
                      pred.z = pred.z,
                      pred.y = as.vector(pred.y))
    }
    
    return(results)
    
  }else if (match.arg(lucid_model) == "parallel"){
    ###LUCID in parallel#####
    if(is.null(colnames(G))){
      Gnames <- paste0("G", 1:ncol(G))
    } else {
      Gnames <- colnames(G)
    }
    colnames(G) <- Gnames


    if(is.null(Z)) {
      stop("Input data 'Z' is missing")
    }
    if(!is.list(Z)) {
      stop("Input data 'Z' should be a list for LUCID in Parallel!")
    } else {
      for(i in 1:length(Z)) {
        if(!is.matrix(Z[[i]])) {
          Z[[i]] <- as.matrix(Z[[i]])
          if(!is.numeric(Z[[i]])) {
            stop("Input data 'Z' should be numeric")
          }
        }
      }}
    Znames <- vector("list", length(Z))
    for(i in 1:length(Z)) {
      if(is.null(colnames(Z[[i]]))){
        Znames[[i]] <- paste0("Z_", i, "_", 1:ncol(Z[[i]]))
      } else {
        Znames[[i]] <- colnames(Z[[i]])
      }}



    N <- nrow(G)
    K <- model$K
    nOmics <- length(K)
    nG <- ncol(G)
    nZ <- as.integer(sapply(Z, ncol))

    dimCoG <- ifelse(is.null(CoG), 0, ncol(CoG))
    dimCoY <- ifelse(is.null(CoY), 0, ncol(CoY))
    # combine G and CoG to adjust for CoG
    if(!is.null(CoG)) {
      G <- cbind(G, CoG)
      Gnames <- c(Gnames, CoGnames)
    }

    modelNames =  model$modelName

    #Get na_pattern, but for prediction, Z should be non-missing
    na_pattern <- vector("list", nOmics)
    for(i in 1:nOmics) {
      na_pattern[[i]] <- check_na(Z[[i]])
      }

    Beta = model$res_Beta$Beta
    Mu = model$res_Mu
    Sigma = model$res_Sigma
    Gamma <- model$res_Gamma$Gamma
    useY_flag <- ifelse(is.null(Y), FALSE, TRUE)
    
    if(g_computation == FALSE){
      Estep_array <- Estep(G = G, Z = Z, Y = Y,
                           Beta = Beta, Mu = Mu, Sigma = Sigma, Delta = Gamma,
                           family = model$family, useY = useY_flag, na_pattern = na_pattern)

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
        pred_X[[i]] <- sapply(1:N, function(x) return(nnet::which.is.max(post.p[[i]][x, ])))
        pred_X[[i]] <- pred_X[[i]] - 1
      }

    # 2 - prediction for Y
    r = Estep_r

    # if 2 omics layers
    if(nOmics == 2) {
      r_matrix <- t(sapply(1:N, function(j) {
        marginSums(lastInd(r,j), margin = i)
      }))
      r_fit <- as.data.frame(r_matrix[, -c(1, K[1] + 1)])
      #Here ip are used to get the weighted Y, we remove reference ip columns as if it's all categorical predicors.
      if(model$family == "gaussian") {
        if(is.null(CoY)) {
          pred_Y <- predict(model$res_Gamma$fit, newdata = r_fit)
        }else{
          r_x <- as.data.frame(cbind(r_fit, CoY))
          colnames(r_x) <- c(paste0("LC", 1:ncol(r_fit)), colnames(CoY))
          pred_Y <- predict(model$res_Gamma$fit, newdata = r_x)

        }}

      if(model$family == "binomial") {
        if(is.null(CoY)) {
          if (response == TRUE){
            pred_Y_prob <- predict(model$res_Gamma$fit, newdata = r_fit, type = "response")
            pred_Y <- ifelse(pred_Y_prob > 0.5, 1, 0)
          }else{pred_Y <- predict(model$res_Gamma$fit, newdata = r_fit, type = "response")}
        }else{
          r_x <- as.data.frame(cbind(r_fit, CoY))
          colnames(r_x) <- c(paste0("LC", 1:ncol(r_fit)), colnames(CoY))
          if (response == TRUE){
            pred_Y_prob <- predict(model$res_Gamma$fit, newdata = r_x, type = "response")
            pred_Y <- ifelse(pred_Y_prob > 0.5, 1, 0)
          }else{pred_Y <- predict(model$res_Gamma$fit, newdata = r_x, type = "response")}
        }
      }}


    # if 3 omics layers
    if(nOmics == 3) {
      r_matrix <- t(sapply(1:N, function(i) {
        c(marginSums(lastInd(r,i), margin = 1),
          marginSums(lastInd(r,i), margin = 2),
          marginSums(lastInd(r,i), margin = 3))
      }))
      r_fit <- as.data.frame(r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1)])

      if(model$family == "gaussian") {
        if(is.null(CoY)) {
          pred_Y <- predict(model$res_Gamma$fit, newdata = r_fit)
        }else{
          r_x <- as.data.frame(cbind(r_fit, CoY))
          colnames(r_x) <- c(paste0("LC", 1:ncol(r_fit)), colnames(CoY))
          pred_Y <- predict(model$res_Gamma$fit, newdata = r_x)

        }}

      if(model$family == "binomial") {
        if(is.null(CoY)) {
          if (response == TRUE){
          pred_Y_prob <- predict(model$res_Gamma$fit, newdata = r_fit, type = "response")
          pred_Y <- ifelse(pred_Y_prob > 0.5, 1, 0)
          }else{pred_Y <- predict(model$res_Gamma$fit, newdata = r_fit, type = "response")}
        }else{
          r_x <- as.data.frame(cbind(r_fit, CoY))
          colnames(r_x) <- c(paste0("LC", 1:ncol(r_fit)), colnames(CoY))
          if (response == TRUE){
          pred_Y_prob <- predict(model$res_Gamma$fit, newdata = r_x, type = "response")
          pred_Y <- ifelse(pred_Y_prob > 0.5, 1, 0)
          }else{pred_Y <- predict(model$res_Gamma$fit, newdata = r_x, type = "response")}
        }
      }}

    # if 4 omics layers
    if(nOmics == 4) {
      r_matrix <- t(sapply(1:N, function(i) {
        c(marginSums(lastInd(r,i), margin = 1),
          marginSums(lastInd(r,i), margin = 2),
          marginSums(lastInd(r,i), margin = 3),
          marginSums(lastInd(r,i), margin = 4))
      }))
      r_fit <- as.data.frame(r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1)])

      if(model$family == "gaussian") {
        if(is.null(CoY)) {
          pred_Y <- predict(model$res_Gamma$fit, newdata = r_fit)
        }else{
          r_x <- as.data.frame(cbind(r_fit, CoY))
          colnames(r_x) <- c(paste0("LC", 1:ncol(r_fit)), colnames(CoY))
          pred_Y <- predict(model$res_Gamma$fit, newdata = r_x)

        }}

      if(model$family == "binomial") {
        if(is.null(CoY)) {
          if (response == TRUE){
            pred_Y_prob <- predict(model$res_Gamma$fit, newdata = r_fit, type = "response")
            pred_Y <- ifelse(pred_Y_prob > 0.5, 1, 0)
          }else{pred_Y <- predict(model$res_Gamma$fit, newdata = r_fit, type = "response")}
        }else{
          r_x <- as.data.frame(cbind(r_fit, CoY))
          colnames(r_x) <- c(paste0("LC", 1:ncol(r_fit)), colnames(CoY))
          if (response == TRUE){
            pred_Y_prob <- predict(model$res_Gamma$fit, newdata = r_x, type = "response")
            pred_Y <- ifelse(pred_Y_prob > 0.5, 1, 0)
          }else{pred_Y <- predict(model$res_Gamma$fit, newdata = r_x, type = "response")}
        }
      }}

    # if 5 omics layers
    if(nOmics == 5) {
      r_matrix <- t(sapply(1:N, function(i) {
        c(marginSums(lastInd(r,i), margin = 1),
          marginSums(lastInd(r,i), margin = 2),
          marginSums(lastInd(r,i), margin = 3),
          marginSums(lastInd(r,i), margin = 4),
          marginSums(lastInd(r,i), margin = 5))
      }))
      r_fit <- as.data.frame(r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1,
                             K[1] + K[2] + K[3] + K[4] + 1)])

      if(model$family == "gaussian") {
        if(is.null(CoY)) {
          pred_Y <- predict(model$res_Gamma$fit, newdata = r_fit)
        }else{
          r_x <- as.data.frame(cbind(r_fit, CoY))
          colnames(r_x) <- c(paste0("LC", 1:ncol(r_fit)), colnames(CoY))
          pred_Y <- predict(model$res_Gamma$fit, newdata = r_x)

        }}

      if(model$family == "binomial") {
        if(is.null(CoY)) {
          if (response == TRUE){
            pred_Y_prob <- predict(model$res_Gamma$fit, newdata = r_fit, type = "response")
            pred_Y <- ifelse(pred_Y_prob > 0.5, 1, 0)
          }else{pred_Y <- predict(model$res_Gamma$fit, newdata = r_fit, type = "response")}
        }else{
          r_x <- as.data.frame(cbind(r_fit, CoY))
          colnames(r_x) <- c(paste0("LC", 1:ncol(r_fit)), colnames(CoY))
          if (response == TRUE){
            pred_Y_prob <- predict(model$res_Gamma$fit, newdata = r_x, type = "response")
            pred_Y <- ifelse(pred_Y_prob > 0.5, 1, 0)
          }else{pred_Y <- predict(model$res_Gamma$fit, newdata = r_x, type = "response")}
        }
      }}
    
    if (g_computation == FALSE){
      results <- list(inclusion.p = post.p,
                      pred.x = pred_X,
                      pred.y = pred_Y)
    }else{
      #compute the weighted (by pred.x and estimated mu from model) predicted mu for each obs of each feature
      pred.z <- vector(mode = "list", length = nOmics)
      for (i in 1:nOmics){
        mu_i = t(Mu[[i]])
        z = ncol(mu_i)
        pred_layer_z = matrix(NA, nrow = N, ncol = z)
        colnames(pred_layer_z) = names(mu_i)
        for (j in 1:N){
          p_j = post.p[[i]][j,]
          mu_j = colSums(p_j * mu_i)
          pred_layer_z[j,] = mu_j
        }
        pred.z[[i]] <- pred_layer_z
      }
      results <- list(inclusion.p = post.p,
                      pred.x = pred_X,
                      pred.z = pred.z,
                      pred.y = pred_Y)
    }
      return(results)
    }
  }







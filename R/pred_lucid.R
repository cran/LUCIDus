############workhorse function for predict_lucid############
############prediction for "early", "parallel"############

pred_lucid <- function(model,
                          lucid_model = c("early", "parallel"),
                          G,
                          Z,
                          Y = NULL,
                          CoG = NULL,
                          CoY = NULL,
                          response = TRUE){
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

    # normalize the log-likelihood to probability
    res.r <- t(apply(res, 1, lse_vec))
    # predicted latent cluster
    pred.x <- sapply(1:n, function(x) return(nnet::which.is.max(res.r[x, ])))
    mu_Y <- cbind(res.r, CoY) %*% gamma$beta


    if(model$family == "normal"){
      pred.y <- mu_Y
    }
    if(model$family == "binary"){
      pred.y <- exp(mu_Y) / (1 + exp(mu_Y))
      if(response == TRUE){
        pred.y <- as.numeric(pred.y > 0.5)
      }
    }

    return(list(inclusion.p = res.r,
                pred.x = pred.x,
                pred.y = as.vector(pred.y)))
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

    Estep_array <- Estep(G = G, Z = Z, Y = Y,
                         Beta = Beta, Mu = Mu, Sigma = Sigma, Delta = Gamma,
                         family = model$family, useY = useY_flag, na_pattern = na_pattern)
    Estep_r <- Estep_to_r(Estep_array = Estep_array,
                          K = K,
                          N = N)

    post.p <- vector(mode = "list", length = nOmics)
    for(i in 1:nOmics) {
      post.p[[i]] = compute_res_r(r = Estep_r, N = N,layer = i)
    }

    # initialize container for predicted value
    pred_X <- vector(mode = "list", length = nOmics)
    pred_z <- vector(mode = "list", length = nOmics)
    # prediction of X and Y based on fitted data


      # 1 - prediction for X
      for (i in 1:nOmics) {
        r_margin <- t(sapply(1:N, function(j) {
          marginSums(lastInd(Estep_r,j), margin = i)
        }))
        pred_X[[i]] <- map(r_margin)
        pred_z[[i]] <- r_margin
      }

    # 2 - prediction for Y
    r = Estep_r

    # if 2 omics layers
    if(nOmics == 2) {
      r_matrix <- t(sapply(1:N, function(j) {
        marginSums(lastInd(r,j), margin = i)
      }))
      r_fit <- as.data.frame(r_matrix[, -c(1, K[1] + 1)])

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

      return(list(inclusion.p = post.p,
                  pred.x = pred_X,
                  pred.y = pred_Y))
    }
  }







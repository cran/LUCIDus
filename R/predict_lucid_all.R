#' @title Predict cluster assignment and outcome based on LUCID model
#'
#' @param model A model fitted and returned by \code{\link{estimate_lucid}}
#' @param lucid_model Specifying LUCID model, "early" for early integration, "parallel" for lucid in parallel
#' "serial" for lucid in serial.
#' @param G Exposures, a numeric vector, matrix, or data frame. Categorical variable
#' should be transformed into dummy variables. If a matrix or data frame, rows
#' represent observations and columns correspond to variables.
#' @param Z Omics data, if "early", an N by M matrix; If "parallel", a list, each element i is a matrix with N rows and P_i features;
#' If "serial", a list, each element i is a matrix with N rows and p_i features or a list with two or more matrices with N rows and a certain number of features
#' @param Y Outcome, a numeric vector. Categorical variable is not allowed. Binary
#' outcome should be coded as 0 and 1.
#' @param CoG Optional, covariates to be adjusted for estimating the latent cluster.
#' A numeric vector, matrix or data frame. Categorical variable should be transformed
#' into dummy variables.
#' @param CoY Optional, covariates to be adjusted for estimating the association
#' between latent cluster and the outcome. A numeric vector, matrix or data frame.
#' Categorical variable should be transformed into dummy variables.
#' @param response If TRUE, when predicting binary outcome, the response will be
#' returned. If FALSE, the linear predictor is returned.
#' @param verbose A flag indicates whether detailed information
#' is printed in console. Default is FALSE.
#' @return A list containing the following components:
#' 1. inclusion.p: A list of inclusion probabilities for each sub-model in the LUCID model.
#' 2. pred.x: A list of predicted values for the data matrix G.
#' 3. pred.y: Predicted values for the response variable Y (if response is TRUE).

#' @export
#'
#' @examples
#' # prepare data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#'
#' # fit lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", K = 2, family = "normal")
#'
#' # prediction on training set
#' pred1 <- predict_lucid(model = fit1, G = G, Z = Z, Y = Y_normal, lucid_model = "early")
#' pred2 <- predict_lucid(model = fit1, G = G, Z = Z, lucid_model = "early")
#' 


predict_lucid <- function(model,
                          lucid_model = c("early", "parallel","serial"),
                          G,
                          Z,
                          Y = NULL,
                          CoG = NULL,
                          CoY = NULL,
                          response = TRUE,
                          verbose = FALSE){

  if (match.arg(lucid_model) == "early" | match.arg(lucid_model) == "parallel"){
    # ========================== Early Integration ==========================
    # ========================== LUCID IN PARALLEL ==========================
    res_pred = pred_lucid(model = model, lucid_model = lucid_model, G = G, Z = Z, Y = Y,
                          CoG = CoG, CoY = CoY, response = response)
    return(res_pred)
  }else if (match.arg(lucid_model) == "serial"){
    # ========================== LUCID IN Serial ==========================
    n <- nrow(G)
    K <- model$K

    ## check data format ==== special for Z  under serial
    if(length(Z) != length(K)) {
      stop("Z and K should be two lists of the same length for LUCID in Serial!")
    }

    if(is.null(Z)) {
      stop("Input data 'Z' is missing")
    }
    if(!is.list(Z)) {
      stop("Input data 'Z' should be a list for LUCID in Serial!")
    }
    else {
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

    # initiate the empty lists to store the predictions for each sub model
    post.p.list <- vector(mode = "list", length = length (K))
    pred.x.list <- vector(mode = "list", length = length (K))

    #loop through each K
    for (i in 1:length(K)){
      if(verbose){
        cat("Predicting LUCID in Serial model",
            paste0("(", "Sub Model Number = ", i, ")"),
            "\n")
      }
      ##Scenario 1: the first serial sub model
      if (i == 1){
        if (is.numeric(K[[1]])){
          #if the first serial sub model is early integration (1 layer)
          temp_pred = pred_lucid(model = model$submodel[[1]], lucid_model = "early", G = G, Z = Z[[1]], Y = NULL,
                                CoG = CoG, CoY = NULL, response = FALSE)

          post.p.list[[1]] = temp_pred$inclusion.p
          pred.x.list[[1]] = temp_pred$pred.x

          post.p = temp_pred$inclusion.p[,-1]

        }else{
          #if the first serial sub model is lucid in parallel
          temp_pred = pred_lucid(model = model$submodel[[1]], lucid_model = "parallel", G = G, Z = Z[[1]], Y = NULL,
                                 CoG = CoG, CoY = NULL, response = FALSE)

          post.p.list[[1]] = temp_pred$inclusion.p
          pred.x.list[[1]] = temp_pred$pred.x

          temp.p = temp_pred$inclusion.p
          temp.p.list = vector(mode = "list", length = length(temp.p))
          for (i in 1:length(temp.p)){
            temp.p.list[[i]] = temp.p[[i]][,-1]
          }
          post.p = matrix(unlist(temp.p.list), nrow = nrow(G), byrow = FALSE)

        }
      }else if (i < length(K)){
        ##Scenario 2: the middle serial sub models
        if (is.numeric(K[[i]])){
          #if the middle serial sub model is early integration (1 layer)
          temp_pred = pred_lucid(model = model$submodel[[i]], lucid_model = "early", G = post.p, Z = Z[[i]], Y = NULL,
                                 CoG = NULL, CoY = NULL, response = FALSE)

          post.p.list[[i]] = temp_pred$inclusion.p
          pred.x.list[[i]] = temp_pred$pred.x

          post.p = temp_pred$inclusion.p[,-1]

        }else{
          #if the first serial sub model is lucid in parallel
          temp_pred = pred_lucid(model = model$submodel[[i]], lucid_model = "parallel", G = post.p, Z = Z[[i]], Y = NULL,
                                 CoG = NULL, CoY = NULL, response = FALSE)

          post.p.list[[i]] = temp_pred$inclusion.p
          pred.x.list[[i]] = temp_pred$pred.x

          temp.p = temp_pred$inclusion.p
          temp.p.list = vector(mode = "list", length = length(temp.p))
          for (i in 1:length(temp.p)){
            temp.p.list[[i]] = temp.p[[i]][,-1]
          }
          post.p = matrix(unlist(temp.p.list), nrow = nrow(G), byrow = FALSE)

        }
      }else if (i == length(K)){
        ##Scenario 3: the last sub model
        if (is.numeric(K[[i]])){
          #if the last serial sub model is early integration (1 layer)
          temp_pred = pred_lucid(model = model$submodel[[i]], lucid_model = "early", G = post.p, Z = Z[[i]], Y = Y,
                                 CoG = NULL, CoY = CoY, response = response)

          post.p.list[[i]] = temp_pred$inclusion.p
          pred.x.list[[i]] = temp_pred$pred.x
          pred.y = temp_pred$pred.y

        }else{
          #if the last serial sub model is parallel (multiple layers)
          temp_pred = pred_lucid(model = model$submodel[[i]], lucid_model = "parallel", G = post.p, Z = Z[[i]], Y = Y,
                                 CoG = NULL, CoY = CoY, response = response)

          post.p.list[[i]] = temp_pred$inclusion.p
          pred.x.list[[i]] = temp_pred$pred.x
          pred.y = temp_pred$pred.y

          }
        }
    }

    return(list(inclusion.p = post.p.list,
                pred.x = pred.x.list,
                pred.y = pred.y))

    }
  }

#' @title Predict Cluster Assignment and Outcome From a Fitted LUCID Model
#' @description Predict cluster assignment and outcome using new data on G, Z, and optional Y.
#' If \code{g_computation = TRUE}, prediction uses only the G-to-X path from the
#' fitted model and returns counterfactual-style predictions under modified G.
#' This function can also be used to extract latent cluster assignments when using
#' the training data as input.
#' @param model A model fitted and returned by \code{\link{estimate_lucid}}
#' @param lucid_model Specifying LUCID model, "early" for early integration,
#' "parallel" for LUCID in parallel, and "serial" for LUCID in serial.
#' @param G Exposures, a numeric vector, matrix, or data frame. Categorical variable
#' should be transformed into dummy variables. If a matrix or data frame, rows
#' represent observations and columns correspond to variables.
#' @param Z Omics data. If "early", an N by M matrix. If "parallel", a list, each
#' element i is a matrix with N rows and P_i features. If "serial", a list, each
#' element i is a matrix with N rows and p_i features (or a list with two or more
#' matrices with N rows and a certain number of features). For
#' \code{g_computation = TRUE}, \code{Z} can be \code{NULL} for "early",
#' "parallel", and "serial".
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
#' @param g_computation If \code{TRUE}, prediction uses only information on G.
#' For this mode, supplied \code{Z} and \code{Y} are ignored for "early",
#' "parallel", and "serial".
#' @param verbose A flag indicates whether detailed information
#' is printed in console. Default is FALSE.
#' @return A list containing:
#' \item{inclusion.p}{Posterior inclusion probabilities for latent clusters (a
#' matrix for "early"; a list by layer for "parallel" and "serial").}
#' \item{pred.x}{Predicted latent-cluster labels (a numeric vector for "early";
#' a list by layer for "parallel" and "serial").}
#' \item{pred.y}{Predicted outcome values. For binary outcomes, this is class
#' labels when \code{response = TRUE} and probabilities when
#' \code{response = FALSE}.}
#' \item{pred.z}{Predicted omics means under g-computation mode
#' (\code{g_computation = TRUE}).}
#' 
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
#' # g-computation style prediction using only G
#' pred_g <- predict_lucid(model = fit1, G = G, Z = NULL, g_computation = TRUE, lucid_model = "early")
#' 


predict_lucid <- function(model,
                          lucid_model = c("early", "parallel","serial"),
                          G,
                          Z = NULL,
                          Y = NULL,
                          CoG = NULL,
                          CoY = NULL,
                          response = TRUE,
                          g_computation = FALSE,
                          verbose = FALSE){
  
  if (g_computation == TRUE){
    if (!is.null(Z) || !is.null(Y)){
      cat("G-computation only uses input for G, and the G-to-X association, input of Z and Y will not be used for prediction.")
    }
  }
  
  if (match.arg(lucid_model) == "early" | match.arg(lucid_model) == "parallel"){
    # ========================== Early Integration ==========================
    # ========================== LUCID IN PARALLEL ==========================
    res_pred = pred_lucid(model = model, lucid_model = lucid_model, G = G, Z = Z, Y = Y,
                          CoG = CoG, CoY = CoY, response = response, g_computation = g_computation,
                          verbose = verbose)
    return(res_pred)
  }else if (match.arg(lucid_model) == "serial"){
    # ========================== LUCID IN Serial ==========================
    n <- nrow(G)
    K <- model$K

    ## check data format ==== special for Z under serial
    # For g-computation, Z is optional and ignored.
    if (!g_computation) {
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
    }

    # initiate the empty lists to store the predictions for each sub model
    post.p.list <- vector(mode = "list", length = length (K))
    pred.x.list <- vector(mode = "list", length = length (K))
    pred.z.list <- vector(mode = "list", length = length (K))
    
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

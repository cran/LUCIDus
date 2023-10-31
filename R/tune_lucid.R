#' @title A wrapper function to perform model selection for LUCID
#'
#' @description Given a grid of K and L1 penalties (incluing Rho_G, Rho_Z_mu and
#' Rho_Z_Cov; for LUCID early only), fit LUCID model over all combinations of K and L1 penalties to
#' determine the optimal penalty. Note that the input of the grid of K differs for different
#' LUCID models. i.e. For LUCID Early, K = 3:5; for LUCID in parallel, K = list(2:3, 2:3);
#' for LUCID in serial, K = list(list(2:3,2),2:3)
#'
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
#' @param K Number of latent clusters. If "early", an integer;If "parallel",a list, each element is an integer/integer vector, same length as Z;
#' If "serial", a list, each element is either an integer or an list of integers, same length as Z.
#' If K is given as a grid, the input of the grid of K differs for different
#' LUCID models. i.e. For LUCID Early, K = 3:5; for LUCID in parallel, K = list(2:3, 2:3);
#' for LUCID in serial, K = list(list(2:3,2),2:3)
#' @param lucid_model Specifying LUCID model, "early" for early integration, "parallel" for lucid in parallel,
#' "serial" for lucid in serial
#' @param family Distribution of outcome. For continuous outcome, use "normal";
#' for binary outcome, use "binary". Default is "normal".
#' @param Rho_G A scalar or a vector. This parameter is the LASSO penalty to regularize
#' exposures. If it is a vector, \code{tune_lucid} will conduct
#' model selection and variable selection. User can try penalties from 0 to 1. Work for LUCID early only.
#' @param Rho_Z_Mu A scalar or a vector. This parameter is the LASSO penalty to
#' regularize cluster-specific means for omics data (Z). If it is a vector,
#' \code{tune_lucid} will conduct model selection and
#' variable selection. User can try penalties from 1 to 100. Work for LUCID early only.
#' @param Rho_Z_Cov A scalar or a vector. This parameter is the graphical LASSO
#' penalty to estimate sparse cluster-specific variance-covariance matrices for omics
#' data (Z). If it is a vector, \code{tune_lucid} will conduct
#' model selection and variable selection. User can try penalties from 0 to 1. Work for LUCID early only.
#' @param verbose_tune A flag to print details of tuning process.
#' @param ... Other parameters passed to \code{estimate_lucid}
#'
#' @export
#'
#' @return A list:
#' \item{best_model}{the best model over different combination of tuning parameters}
#' \item{tune_list}{a data frame contains combination of tuning parameters and c
#' orresponding BIC}
#' \item{res_model}{a list of LUCID models corresponding to each combination of
#' tuning parameters}
#'
#' @examples
#' \dontrun{
#' # use simulated data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#'
#' # find the optimal model over the grid of K
#' tune_K <- tune_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early",
#'  useY = FALSE, tol = 1e-3,
#' seed = 1, K = 2:5)
#'
#' # tune penalties
#' tune_Rho_G <- tune_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early",
#'  useY = FALSE, tol = 1e-3,
#' seed = 1, K = 2, Rho_G = c(0.1, 0.2, 0.3, 0.4))
#' tune_Rho_Z_Mu <- tune_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", 
#' useY = FALSE, tol = 1e-3,
#' seed = 1, K = 2, Rho_Z_Mu = c(10, 20, 30, 40))
#' tune_Rho_Z_Cov <- tune_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", 
#' useY = FALSE, tol = 1e-3,
#' seed = 1, K = 2, Rho_Z_Cov = c(0.1, 0.2, 0.3))
#'}
#'
tune_lucid <- function(G,
                       Z,
                       Y,
                       CoG = NULL,
                       CoY = NULL,
                       family = c("normal","binary"),
                       K,
                       lucid_model = c("early", "parallel","serial"),
                       Rho_G = 0,
                       Rho_Z_Mu = 0,
                       Rho_Z_Cov = 0,
                       verbose_tune = FALSE,
                       ...){
  if (match.arg(lucid_model) == "early" | match.arg(lucid_model) == "parallel"){
      tune_lucid_auxi(G = G,
                      Z = Z,
                      Y = Y,
                      CoG = CoG,
                      CoY = CoY,
                      family = family,
                      K = K,
                      lucid_model = lucid_model,
                      Rho_G = Rho_G,
                      Rho_Z_Mu = Rho_Z_Mu,
                      Rho_Z_Cov = Rho_Z_Cov,
                      verbose_tune = verbose_tune,
                      ...)
  }else if (match.arg(lucid_model) == "serial"){
    #tune lucid for parallel, Rho disabled, you can only have a fixed Rho for the LUCID early sub models
    if (length(Rho_G) > 1 | length(Rho_Z_Mu) > 1 | length(Rho_Z_Cov) > 1){
      stop("Tune LUCID in serial can't tune for multiple regualrities for now!")
    }
    #construct the tune list given K
    ##check if all the elements in K_list[[i]] are not lists, if so, it's a serial model of only early sub model
    if (all(!sapply(K, is.list))){
      K_list = expand_list_grid(K)
      for (i in 1:length(K_list)){
        K_list[[i]] <- K_list[[i]][-length(K_list[[i]])]
        all_not_lists <- all(!sapply(K_list[[i]], is.list))

        if (all_not_lists) {
        K_list[[i]] = unlist(K_list[[i]])
      }
    }}else{
      #There are parallel sub models
      K_list = expand_list_grid(K)
      for (i in 1:length(K_list)){
        K_list[[i]] <- K_list[[i]][-length(K_list[[i]])]

        for (j in 1:length(K_list[[i]])){
          if (length(K_list[[i]][[j]])>1){
            names(K_list[[i]][[j]]) = NULL
            K_list[[i]][[j]] <- as.list(K_list[[i]][[j]])
          }
        }
      }}


    K_matrix = as.matrix(K_list)


    # use bic to conduct model selection
    bic <- rep(0, length(K_list))
    model_list <- vector(mode = "list", length = length(K_list))
    for(i in 1:length(K_list)) {
      model_list[[i]] <- estimate_lucid(G = G, Z = Z, Y = Y,
                                        CoG = CoG,
                                        CoY = CoY,
                                        family = family,
                                        lucid_model = "serial",
                                        K = K_list[[i]],
                                        Rho_G = Rho_G,
                                        Rho_Z_Mu = Rho_Z_Mu,
                                        Rho_Z_Cov = Rho_Z_Cov,
                                        verbose = verbose_tune,
                                        ...)
      bic[i] <- cal_bic_serial(model_list[[i]])
    }
    model_opt_index <- which(bic == min(bic))
    K_matrix <- cbind(K_matrix, bic)
    return(list(tune_K = K_matrix,
                model_list = model_list,
                model_opt = model_list[[model_opt_index]]))

  }

  }



tune_lucid_auxi <- function(G,
                       Z,
                       Y,
                       CoG = NULL,
                       CoY = NULL,
                       family = "normal",
                       K = 2:5,
                       lucid_model = c("early", "parallel"),
                       Rho_G = 0,
                       Rho_Z_Mu = 0,
                       Rho_Z_Cov = 0,
                       verbose_tune = FALSE,
                       ...){

  if (match.arg(lucid_model) == "early"){
  #tune lucid for early, Rho enabled
  # combinations of tuning parameters
  tune_list <- expand.grid(K, Rho_G, Rho_Z_Mu, Rho_Z_Cov)
  colnames(tune_list) <- c("K", "Rho_G", "Rho_Z_Mu", "Rho_Z_Cov")
  m <- nrow(tune_list)
  tune_list$BIC <- rep(0, m)
  if(m > 1) {
    cat("Tuning LUCID model \n \n")
  }

  # fit models for each combination
  res_model <- vector(mode = "list",
                      length = m)
  for(i in 1:m) {
    fit <- try(estimate_lucid(G = G,
                         Z = Z,
                         Y = Y,
                         CoG = CoG,
                         CoY = CoY,
                         family = family,
                         lucid_model = "early",
                         K = tune_list[i, 1],
                         Rho_G = tune_list[i, 2],
                         Rho_Z_Mu = tune_list[i, 3],
                         Rho_Z_Cov = tune_list[i, 4],
                         verbose = verbose_tune,
                         ...))
    if("try-error" %in% class(fit)) {
      tune_list[i, 5] <- NA
    } else {
      tune_list[i, 5] <- summary_lucid(fit)$BIC
    }
    res_model[[i]] <- fit
  }
  x <- min(tune_list[, 5], na.rm = TRUE)
  if(is.na(x)) {
    stop("LUCID model fails to converge given current tuning parameters")
  }
  best_model <- res_model[[which(tune_list[, 5]== x)]]
  return(list(best_model = best_model,
              tune_list = tune_list,
              res_model = res_model))
  }

######LUCID in parallel, needs work (from model_select)....######
  else if (match.arg(lucid_model) == "parallel"){
      #tune lucid for parallel, Rho disabled
      # use bic to conduct model selection
      if (length(Rho_G) > 1 | length(Rho_Z_Mu) > 1 | length(Rho_Z_Cov) > 1){
        stop("Tune LUCID in parallel can't tune for multiple regualrities for now!")
      }
      K_list = K
      # change K_list into a matrix
      K_matrix <- as.matrix(expand.grid(K_list))
      if(min(K_matrix) < 2) {
        stop("minimum K should be 2")
      }

      bic <- rep(0, nrow(K_matrix))
      model_list <- vector(mode = "list", length = nrow(K_matrix))
      for(i in 1:nrow(K_matrix)) {
        model_list[[i]] <- estimate_lucid(G = G, Z = Z, Y = Y,
                                          CoG = CoG,
                                          CoY = CoY,
                                          family = family,
                                          lucid_model = "parallel",
                                          K = K_matrix[i, ],
                                          verbose = verbose_tune,
                                          ...)
        bic[i] <- cal_bic_parallel(model_list[[i]])
      }
      model_opt_index <- which(bic == min(bic))
      K_matrix <- cbind(K_matrix, bic)
      return(list(tune_K = K_matrix,
                  model_list = model_list,
                  model_opt = model_list[[model_opt_index]]))

  }
}

#expand list of list of integers or integers for lucid in serial
expand_list_grid <- function(lst) {
  if (length(lst) == 0) {
    return(list(list(list())))
  }

  # Extract the first element from the list
  first_element <- lst[[1]]

  # Remove the first element from the list
  rest_of_list <- lst[-1]

  # Recursively expand the rest of the list
  expanded_rest <- expand_list_grid(rest_of_list)

  # Create a grid using expand.grid for the first element
  grid <- expand.grid(first_element)

  # Combine the grid with each element of the expanded rest of the list
  result <- list()
  for (i in 1:nrow(grid)) {
    for (j in 1:length(expanded_rest)) {
      combined_element <- c(list(unlist(grid[i, ])), expanded_rest[[j]])
      result <- c(result, list(combined_element))
    }
  }

  return(result)
}


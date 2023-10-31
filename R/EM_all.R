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
#' @param K Number of latent clusters. If "early", an integer greater or equal to 2; If "parallel",an integer vector, same length as Z, with each element being an interger greater or equal to 2;
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
#' exposures. If user wants to tune the penalty, use the wrapper
#' function \code{lucid}. Now only achieved for LUCID early integration.
#' @param Rho_Z_Mu A scalar. This parameter is the LASSO penalty to
#' regularize cluster-specific means for omics data (Z). If user wants to tune the
#' penalty, use the wrapper function \code{lucid}.Now only achieved for LUCID early integration.
#' @param Rho_Z_Cov A scalar. This parameter is the graphical LASSO
#' penalty to estimate sparse cluster-specific variance-covariance matrices for omics
#' data (Z). If user wants to tune the penalty, use the wrapper function \code{lucid}.
#' Now only achieved for LUCID early integration.
#' @param family The distribution of the outcome
#' @param seed Random seed to initialize the EM algorithm
#' @param init_impute Method to initialize the imputation of missing values in
#' LUCID. \code{mix} will use \code{mclust:imputeData} to implement EM Algorithm
#' for Unrestricted General Location Model by the mix package to impute the missing values in omics
#' data; \code{lod} will initialize the imputation via replacing missing values by
#' LOD / sqrt(2). LOD is determined by the minimum of each variable in omics data.
#' @param init_par For "early", an interface to initialize EM algorithm, if mclust,
#' initiate the parameters using the \code{mclust} package, if random, initiate the parameters
#' by drawing from a uniform distribution;
#' For "parallel", mclust is the default for quick convergence;
#' For "serial", each sub-model follows the above depending on it is a "early" or "parallel"
#' @param verbose A flag indicates whether detailed information for each iteration
#' of EM algorithm is printed in console. Default is FALSE.
#'
#' @import mclust
#' @import stats
#' @import utils
#' @import glasso
#' @import glmnet
#' @return A list contains the object below:
#' 1. res_Beta: estimation for G->X associations
#' 2. res_Mu: estimation for the mu of the X->Z associations
#' 3. res_Sigma: estimation for the sigma of the X->Z associations
#' 4. res_Gamma: estimation for X->Y associations
#' 5. inclusion.p: inclusion probability of cluster assignment for each observation
#' 6. K: umber of latent clusters for "early"/list of numbers of latent clusters for "parallel" and "serial"
#' 7. var.names: names for the G, Z, Y variables
#' 8. init_omic.data.model: pre-specified geometric model of multi-omics data
#' 9. likelihood: converged LUCID model log likelihood
#' 10. family: the distribution of the outcome
#' 11. select: for LUCID early integration only, indicators of whether each exposure and omics feature is selected 
#' 12. useY: whether this LUCID model is supervised
#' 13. Z: multi-omics data
#' 14. init_impute: pre-specified imputation method
#' 15. init_par: pre-specified parameter initialization method
#' 16. Rho: for LUCID early integration only, pre-specified regularity tuning parameter 
#' 17. N: number of observations
#' 18. submodel: for LUCID in serial only, storing all the submodels
#' @examples
#' i <- 1008
#' set.seed(i)
#' G <- matrix(rnorm(500), nrow = 100)
#' Z1 <- matrix(rnorm(1000),nrow = 100)
#' Z2 <- matrix(rnorm(1000), nrow = 100)
#' Z3 <- matrix(rnorm(1000), nrow = 100)
#' Z4 <- matrix(rnorm(1000), nrow = 100)
#' Z5 <- matrix(rnorm(1000), nrow = 100)
#' Z <- list(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4, Z5 = Z5)
#' Y <- rnorm(100)
#' CoY <- matrix(rnorm(200), nrow = 100)
#' CoG <- matrix(rnorm(200), nrow = 100)
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y, K = list(2,2,2,2,2),
#' lucid_model = "serial",
#' family = "normal",
#' seed = i,
#' CoG = CoG, CoY = CoY,
#' useY = TRUE)
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
                      init_impute = c("mix", "lod"),
                      init_par = c("mclust", "random"),
                      verbose = FALSE) {
  if (match.arg(lucid_model) == "early" | match.arg(lucid_model) == "parallel"){
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
                         verbose = verbose)
    return(results)
  }else{
    # ========================== LUCID IN Serial ==========================
    ## 1.1 check data format ==== special for Z  under serial
    #The Z input for LUCID in serial should be a list
    #length(Z) = length (K)
    #for each element of K that is a list
    #the corresponding element of Z should also be a list

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

    # initiate the empty lists to store the parameter estimates
    post.p.list <- vector(mode = "list", length = length (K))
    res.mu.list <- vector(mode = "list", length = length (K))
    res.sigma.list <- vector(mode = "list", length = length (K))
    #delta is for association between LCs for each sub models
    res.delta.list <- vector(mode = "list", length = length (K)-1)
    Znames <- vector(mode = "list", length = length (K))
    submodel <- vector(mode = "list", length = length (K))

    #loop through each K
    for (i in 1:length(K)){
      if(verbose){
        cat("Fitting LUCID in Serial model",
            paste0("(", "Sub Model Number = ", i, ")"),
            "\n")
      }
      set.seed(seed + i * 1900)
      #simulate random Y cause we don't need Y except the last sub model
      Y_rand = runif(nrow(G))
      ##Scenario 1: the first serial sub model
      if (i == 1){
      if (is.numeric(K[[1]])){
       #if the first serial sub model is early integration (1 layer)
        temp_model = est_lucid(lucid_model = "early",
                               G = G,
                               Z = Z[[1]],
                               Y = Y_rand,
                               CoG = CoG,
                               CoY = NULL, K = unlist(K[[1]]),
                               init_omic.data.model = init_omic.data.model,
                               useY = FALSE,
                               tol = tol,
                               max_itr = max_itr,
                               max_tot.itr = max_tot.itr,
                               Rho_G = Rho_G,
                               Rho_Z_Mu = Rho_Z_Mu,
                               Rho_Z_Cov = Rho_Z_Cov,
                               family = "normal",
                               seed = seed,
                               init_impute = init_impute,
                               init_par = init_par,
                               verbose = verbose)
        #update parameters
        post.p.list[[1]] = temp_model$inclusion.p
        res.mu.list[[1]] = temp_model$res_Mu
        res.sigma.list[[1]] = temp_model$res_Sigma
        res_Beta = temp_model$res_Beta
        Gnames = temp_model$var.names$Gnames
        Znames[[1]] = temp_model$var.names$Znames
        submodel[[1]] = temp_model
        #update PIP as the input for the next sub model, excluding the PIP for the reference cluster
        post.p = temp_model$inclusion.p[,-1]

      }else{
        #if the first serial sub model is lucid in parallel
        temp_model = est_lucid(lucid_model = "parallel",
                               G = G,
                               Z = Z[[1]],
                               Y = Y_rand,
                               CoG = CoG,
                               CoY = NULL, K = unlist(K[[1]]),
                               init_omic.data.model = init_omic.data.model,
                               useY = FALSE,
                               tol = tol,
                               max_itr = max_itr,
                               max_tot.itr = max_tot.itr,
                               Rho_G = Rho_G,
                               Rho_Z_Mu = Rho_Z_Mu,
                               Rho_Z_Cov = Rho_Z_Cov,
                               family = "normal",
                               seed = seed,
                               init_impute = init_impute,
                               init_par = init_par,
                               verbose = verbose)
        #update parameters
        post.p.list[[1]] = temp_model$inclusion.p
        res.mu.list[[1]] = temp_model$res_Mu
        res.sigma.list[[1]] = temp_model$res_Sigma
        res_Beta = temp_model$res_Beta
        Gnames = temp_model$var.names$Gnames
        Znames[[1]] = temp_model$var.names$Znames
        submodel[[1]] = temp_model
        #update PIP as the input for the next sub model, excluding the PIP for the reference cluster
        temp.p = temp_model$inclusion.p
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
          temp_model = est_lucid(lucid_model = "early",
                                 G = post.p,
                                 Z = Z[[i]],
                                 Y = Y_rand,
                                 CoG = NULL,
                                 CoY = NULL, K = unlist(K[[i]]),
                                 init_omic.data.model = init_omic.data.model,
                                 useY = FALSE,
                                 tol = tol,
                                 max_itr = max_itr,
                                 max_tot.itr = max_tot.itr,
                                 Rho_G = Rho_G,
                                 Rho_Z_Mu = Rho_Z_Mu,
                                 Rho_Z_Cov = Rho_Z_Cov,
                                 family = "normal",
                                 seed = seed,
                                 init_impute = init_impute,
                                 init_par = init_par,
                                 verbose = verbose)
          #update parameters
          post.p.list[[i]] = temp_model$inclusion.p
          res.mu.list[[i]] = temp_model$res_Mu
          res.sigma.list[[i]] = temp_model$res_Sigma
          res.delta.list[[i-1]] = temp_model$res_Beta
          Znames[[i]] = temp_model$var.names$Znames
          submodel[[i]] = temp_model
          #update PIP as the input for the next sub model, excluding the PIP for the reference cluster
          post.p = temp_model$inclusion.p[,-1]
        }else{
          #if the middle serial sub model is parallel (multiple layers)
          temp_model = est_lucid(lucid_model = "parallel",
                                 G = post.p,
                                 Z = Z[[i]],
                                 Y = Y_rand,
                                 CoG = NULL,
                                 CoY = NULL, K = unlist(K[[i]]),
                                 init_omic.data.model = init_omic.data.model,
                                 useY = FALSE,
                                 tol = tol,
                                 max_itr = max_itr,
                                 max_tot.itr = max_tot.itr,
                                 Rho_G = Rho_G,
                                 Rho_Z_Mu = Rho_Z_Mu,
                                 Rho_Z_Cov = Rho_Z_Cov,
                                 family = "normal",
                                 seed = seed,
                                 init_impute = init_impute,
                                 init_par = init_par,
                                 verbose = verbose)
          #update parameters
          post.p.list[[i]] = temp_model$inclusion.p
          res.mu.list[[i]] = temp_model$res_Mu
          res.sigma.list[[i]] = temp_model$res_Sigma
          res.delta.list[[i-1]] = temp_model$res_Beta
          Znames[[i]] = temp_model$var.names$Znames
          submodel[[i]] = temp_model
          #update PIP as the input for the next sub model, excluding the PIP for the reference cluster
          temp.p = temp_model$inclusion.p
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
          temp_model = est_lucid(lucid_model = "early",
                                 G = post.p,
                                 Z = Z[[i]],
                                 Y = Y,
                                 CoG = NULL,
                                 CoY = CoY, K = unlist(K[[i]]),
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
                                 verbose = verbose)
          #update parameters
          post.p.list[[i]] = temp_model$inclusion.p
          res.mu.list[[i]] = temp_model$res_Mu
          res.sigma.list[[i]] = temp_model$res_Sigma
          res.delta.list[[i-1]] = temp_model$res_Beta
          res_Gamma = temp_model$res_Gamma
          Znames[[i]] = temp_model$var.names$Znames
          Ynames = temp_model$var.names$Ynames
          submodel[[i]] = temp_model

        }else{
          #if the last serial sub model is parallel (multiple layers)
          temp_model = est_lucid(lucid_model = "parallel",
                                  G = post.p,
                                  Z = Z[[i]],
                                  Y = Y,
                                  CoG = NULL,
                                  CoY = CoY, K = unlist(K[[i]]),
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
                                  verbose = verbose)
          #update parameters
          post.p.list[[i]] = temp_model$inclusion.p
          res.mu.list[[i]] = temp_model$res_Mu
          res.sigma.list[[i]] = temp_model$res_Sigma
          res.delta.list[[i-1]] = temp_model$res_Beta
          res_Gamma = temp_model$res_Gamma
          Znames[[i]] = temp_model$var.names$Znames
          Ynames = temp_model$var.names$Ynames
          submodel[[i]] = temp_model
        }
      }
    }
    if(verbose){
      cat("Success: LUCID in Serial Model is constructed!", "\n\n")
    }
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
                    #likelihood = loglik_update, *??? needs to discuss
                    inclusion.p = post.p.list,
                    family = family,
                    #select = list(selectG = selectG, selectZ = selectZ),
                    useY = useY,
                    Z = Z,
                    #z = Estep_r,
                    init_impute = init_impute,
                    init_par = init_par,
                    submodel = submodel
                    #Rho = list(Rho_G = Rho_G,
                    #Rho_Z_Mu = Rho_Z_Mu,
                    #Rho_Z_Cov = Rho_Z_Cov)
    )
    class(results) <- c("lucid_serial")
    return(results)
  }
}

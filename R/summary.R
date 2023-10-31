#' @title Summarize results of the LUCID model
#'
#' @param object A LUCID model fitted by \code{\link{estimate_lucid}}
#' @param boot.se An object returned by \code{\link{boot_lucid}},
#' which contains the bootstrap confidence intervals
#' @return A list, containing the extracted key parameters from the LUCID model that can be used to print the summary table
#' @export
#' @examples
#' \donttest{
#' # use simulated data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#'
#' # fit lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", family = "normal", K = 2,
#' seed = 1008)
#'
#' # conduct bootstrap resampling
#' boot1 <- boot_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", model = fit1, R = 100)
#'
#' # summarize lucid model
#' summary_lucid(fit1)
#'
#' # summarize lucid model with bootstrap CIs
#' summary_lucid(fit1, boot.se = boot1)
#' }

summary_lucid <- function(object, boot.se = NULL){
  if (inherits(object, "early_lucid") | inherits(object, "lucid_parallel")){
    summary_lucid_auxi(object = object, boot.se = boot.se)
  }
  else if (inherits(object, "lucid_serial")){
    K = object$K
    submodels = object$submodel
    n_submodels = length(submodels)
    summary.list <- vector(mode = "list", length = n_submodels)
    for (i in 1:n_submodels){
      summary.list[[i]] = summary_lucid_auxi(submodels[[i]])

      # correct the submodel elements
      #if (i != 1){
        #names(summary.list[[i]])[1] = "delta"
      #}
      #if (i != n_submodels){
        #summary.list[[i]] = summary.list[[i]][-3]
      #}
    }
    BIC = cal_bic_serial(object)
    loglik = cal_loglik_serial(object)
    results <- list(summary.list = summary.list,
                    BIC = BIC,
                    loglik = loglik
                    )

    #return a list of sumlucid object for each submodel
    class(results) <- "sumlucid_serial"
    return(results)
  }
}


summary_lucid_auxi <- function(object, boot.se = NULL){
  if (inherits(object, "early_lucid")){
    s1 <- object$select$selectG
    s2 <- object$select$selectZ
    nG <- sum(s1)
    nZ <- sum(s2)
    K <- object$K
    gamma <- object$res_Gamma$beta
    #obtain number of parameters
    if(object$family == "normal"){
      nY <- length(object$res_Gamma$beta) + length(object$res_Gamma$sigma)
    }
    if(object$family == "binary"){
      nY <- length(object$res_Gamma$beta)
    }
    npars <- (nG + 1) * (K - 1) + (nZ * K + nZ^2 * K) + nY
    BIC <- -2 * object$likelihood + npars * log(nrow(object$inclusion.p))
    results <- list(beta = object$res_Beta[, c(TRUE, s1)],
                    mu = object$res_Mu[, s2],
                    gamma = object$res_Gamma,
                    family = object$family,
                    K = K,
                    BIC = BIC,
                    loglik = object$likelihood,
                    boot.se = boot.se)
    class(results) <- "sumlucid_early"
    return(results)
  }
  if (inherits(object, "lucid_parallel")){
    ##not having regularity yet, to be added
    s1 <- object$select$selectG
    s2 <- object$select$selectZ
    nG <- sum(s1)
    nZ <- sapply(s2,sum)
    K <- object$K
    #obtain number of parameters
    if(object$family == "gaussian"){
      nY <- length(object$res_Gamma$Gamma$mu) + length(object$res_Gamma$Gamma$sd)
    }
    if(object$family == "binomial"){
      #binary summary res_Gamma$Gamma$mu is unclear, use object$res_Gamma$fit$coefficients instead
      nY <- length(object$res_Gamma$fit$coefficients)
    }
    #initiate number of parameters
    npars = 0
    #compute number of parameters for G to X association
    for (i in 1:length(K)){
      npars_new = (nG + 1) * (K[i] - 1)
      npars = npars + npars_new
    }
    #compute number of parameters for X to Z association and add
    for (i in 1:length(K)){
      npars_new = (nZ[i] * K[i] + nZ[i] * nZ[i] * K[i])
      npars = npars + npars_new
    }
    #compute number of parameters for X to Y association and add
    npars = npars + nY

    BIC <- -2 * object$likelihood + npars * log(nrow(object$inclusion.p[[1]]))

    results <- list(beta = object$res_Beta,
                    mu = object$res_Mu,
                    Gamma = object$res_Gamma,
                    family = object$family,
                    K = K,
                    BIC = BIC,
                    loglik = object$likelihood,
                    #BOOT.SE IS NULL FOR NOW
                    boot.se = NULL
    )
    class(results) <- "sumlucid_parallel"
    return(results)

  }
}

#' @title Print the output of LUCID in a nicer table
#'
#' @param x An object returned by \code{summary_lucid}
#' @param ... Other parameters to be passed to \code{print.sumlucid}
#' @return A nice table/several nice tables of the summary of the LUCID model
#' @export print.sumlucid
#' @export 
#' @examples
#' \donttest{
#' # use simulated data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#'
#' # fit lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", family = "normal", K = 2,
#' seed = 1008)
#'
#' # conduct bootstrap resampling
#' boot1 <- boot_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", model = fit1, R = 100)
#'
#' # print the summary of the lucid model in a table
#' print.sumlucid(summary_lucid(fit1))
#'
#' # print the summary of the lucid model with bootstrap CIs in a table
#' print.sumlucid(summary_lucid(fit1, boot.se = boot1))
#' }

print.sumlucid<- function(x, ...){
  if (inherits(x, "sumlucid_early") | inherits(x, "sumlucid_parallel")){
    print.sumlucid_auxi(x)
  }
  else if (inherits(x, "sumlucid_serial")){
    sum_list = x$summary.list
    for (i in 1:length(sum_list)){
      if (i == 1){
        print.auxi.serial.scen1(sum_list[[i]])
      }else if (i != length(sum_list)){
        print.auxi.serial.scen2(sum_list[[i]])
      }else if (i == length(sum_list)){
        print.auxi.serial.scen3(sum_list[[i]])
      }
    }
    cat("\n \n")
    cat("----------Overall Summary of the LUCID in Serial model---------- \n \n")
    cat("log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")


  }
}


print.sumlucid_auxi <- function(x, ...){
  if(inherits(x, "sumlucid_early")){
  K <- x$K
  beta <- as.data.frame(x$beta)
  dim1 <- ncol(beta) - 1
  z.mean <- as.data.frame(t(x$mu))
  cat("----------Summary of the LUCID Early Integration model---------- \n \n")
  cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")
  y <- switch(x$family, normal = f.normal.early,
              binary = f.binary.early)
  y(x$gamma, K, se = x$boot.se$gamma)
  cat("\n")
  cat("(2) Z: mean of omics data for each latent cluster \n")
  if(is.null(x$boot.se)){
    colnames(z.mean) <- paste0("mu_cluster", 1:K)
    print(z.mean)
  } else{
    print(x$boot.se$mu)
  }
  cat("\n")
  cat("(3) E: odds ratio of being assigned to each latent cluster for each exposure \n")
  if(is.null(ncol(beta))) {
    cat("no exposure is selected given current penalty Rho_G, please use a smaller penalty")
  } else {
    dd <- as.matrix(as.data.frame(beta)[2:K, 2:ncol(beta)])
    g.or <- data.frame(beta = unlist(split(dd, row(dd))))
    rownames(g.or) <- paste0(colnames(beta)[-1], ".cluster", sapply(2:K, function(x) return(rep(x, dim1))))
    if(is.null(x$boot.se)){
      g.or$OR <- exp(g.or$beta)
      print(g.or)
    } else{
      print(x$boot.se$beta)
      }
    }
  }
  if(inherits(x, "sumlucid_parallel")){
    K <- x$K
    #list of betas for each layer
    beta <- x$beta$Beta
    dim1 <- ncol(x$beta$Beta[[1]]) - 1
    #list of Z means for each layer, transposed
    z.mean <- x$mu
    cat("----------Summary of the LUCID in Parallel model---------- \n \n")
    cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")
    y <- switch(x$family, gaussian = f.normal.parallel,
                binomial = f.binary.parallel)
    y(x$Gamma, K, se = x$boot.se$gamma)
    cat("\n")
    cat("(2) Z: mean of omics data for each latent cluster of each layer \n")
    if(is.null(x$boot.se)){
      for (i in 1:length(K)){
        cat("Layer ",i, "\n \n")
        colnames(z.mean[[i]]) <- paste0("mu_cluster", 1:K[i])
        print(z.mean[[i]])
        cat("\n \n")
      }

    }else{
      print(x$boot.se$mu)
    }
    cat("\n")
    cat("(3) E: odds ratio of being assigned to each latent cluster for each exposure for each layer \n")
    if(is.null(beta)) {
      cat("no exposure is selected given current penalty Rho_G, please use a smaller penalty")
    } else {
      for (i in 1:length(K)){
      cat("Layer ",i, "\n \n")
      dd <- as.matrix(as.data.frame(beta[[i]][, -1]))
      g.or <- data.frame(beta = unlist(split(dd, row(dd))))

      rownames(g.or) <- paste0(colnames(beta[[i]])[-1], ".cluster", sapply(2:K[i], function(x) return(rep(x, dim1))))


      if(is.null(x$boot.se)){
        g.or$OR <- exp(g.or$beta)
        print(g.or)
      } else{
        print(x$boot.se$beta)
      }
      cat("\n \n")
      }
    }
  }
}


# summarize output of normal outcome for early
f.normal.early <- function(x, K, se){

  cat("(1) Y (continuous outcome): effect size of Y for each latent cluster (and effect of covariates if included) \n")

  if(!is.null(se)){
    y <- se
  } else {
    beta <- x$beta
    y <- as.data.frame(beta)
    #reference
    y[1,] = 0
    row.names(y)[1:K] <- paste0("cluster", 1:K)
    colnames(y) <- "Gamma"

  }
  print(y)
}


# summarize output of binary outcome for early
f.binary.early <- function(x, K, se){
  cat("(1) Y (binary outcome): log odds of Y for cluster 1 (reference) and log OR for rest cluster (and log OR of covariate if included)\n")
  gamma <- as.data.frame(x$beta)
  gamma[1,] = 0
  colnames(gamma) <- "gamma"
  if(is.null(se)){
    gamma$`exp(gamma)` <- exp(gamma$gamma)
  } else{
    gamma <- cbind(gamma, se[, -1])
  }
  print(gamma)
}

# summarize output of normal outcome for parallel
f.normal.parallel <- function(x, K, se){

  cat("(1) Y (continuous outcome): effects of each non-reference latent cluster for each layer of Y (and effect of covariates if included) \n")

  if(!is.null(se)){
    y <- se
  } else {
    gamma <- x$fit$coefficients[-1]
    gamma <- as.data.frame(gamma)
    counter = 0
    for (i in 1:length(K)){
      non_ref_num = K[i] -1
      end = counter + non_ref_num
      counter = counter + 1
      row.names(gamma)[counter:end] <- paste0("cluster", 2:K[i],"Layer",i)
      counter = end
    }

    colnames(gamma) <- "Gamma"
  }
  print(gamma)
}


# summarize output of binary outcome for parallel
f.binary.parallel <- function(x, K, se){
  cat("(1) Y (binary outcome): log OR for non-reference clusters for each layer (and log OR of covariate if included)\n")
  gamma <- as.data.frame(x$fit$coefficients[-1])
  colnames(gamma) <- "gamma"
  counter = 0
  for (i in 1:length(K)){
    non_ref_num = K[i] -1
    end = counter + non_ref_num
    counter = counter + 1
    row.names(gamma)[counter:end] <- paste0("cluster", 2:K[i],"Layer",i)
    counter = end
  }
  if(is.null(se)){
    gamma$`exp(gamma)` <- exp(gamma$gamma)
  } else{
    gamma <- cbind(gamma, se[, -1])
  }
  print(gamma)
}

##########functions for LUCID in parallel##########




# rearrange cluster order
#
# for continuous outcome - use the cluster combination corresponding to smallest
# mean as the reference cluster
get_ref_cluster <- function(Delta) {
  K <- Delta$K
  mu <- Delta$mu
  mu_matrix <- vec_to_array(K = K, mu = mu)
  ref_index <- which(mu_matrix == min(mu_matrix))
  ref <- arrayInd(ref_index, .dim = K)
  return(ref)
}


# re-arrange parameters for Delta
reorder_Delta <- function(ref, Delta) {
  K <- Delta$K
  mu_matrix <- vec_to_array(K = K, mu = Delta$mu)
  mu <- mu_matrix[ref]

  # if 1 omics layers
  if(length(K) == 1) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1] - mu_matrix[ref[1], 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])


    K_order <- list(K1 = k1_order)
  }

  # if 2 omics layers
  if(length(K) == 2) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1] - mu_matrix[ref[1], 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i] - mu_matrix[1, ref[2]]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order)
  }


  # if 3 omics layers
  if(length(K) == 3) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1] - mu_matrix[ref[1], 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1] - mu_matrix[1, ref[2], 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i] - mu_matrix[1, 1, ref[3]]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])


    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order)
  }


  # if 4 omics layers
  if(length(K) == 4) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1,1] - mu_matrix[ref[1], 1, 1,1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1,1] - mu_matrix[1, ref[2], 1,1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i,1] - mu_matrix[1, 1, ref[3],1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])

    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i] - mu_matrix[1, 1, 1, ref[4]]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order)
  }

  # if 5 omics layers
  if(length(K) == 5) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1, 1, 1] - mu_matrix[ref[1], 1, 1, 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1, 1, 1] - mu_matrix[1, ref[2], 1, 1, 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i, 1, 1] - mu_matrix[1, 1, ref[3], 1, 1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])

    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i, 1] - mu_matrix[1, 1, 1, ref[4], 1]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])

    # reorder K5
    mu_k5 <- rep(0, K[5])
    for(i in 1:K[5]) {
      mu_k5[i] <- mu_matrix[1, 1, 1, 1, i] - mu_matrix[1, 1, 1, 1, ref[5]]
    }
    mu_k5_sort <- sort(mu_k5)
    mu <- c(mu, mu_k5_sort[mu_k5_sort != 0])
    # order of re-arranged cluster for omics 5
    k5 <- order(mu_k5)
    k5_order <- c(ref[5], k5[k5 != ref[5]])


    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order,
                    K5 = k5_order)
  }


  # if 6 omics layers
  if(length(K) == 6) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1, 1, 1, 1] - mu_matrix[ref[1], 1, 1, 1, 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])

    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1, 1, 1, 1] - mu_matrix[1, ref[2], 1, 1, 1, 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])

    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i, 1, 1, 1] - mu_matrix[1, 1, ref[3], 1, 1, 1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])

    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i, 1, 1] - mu_matrix[1, 1, 1, ref[4], 1, 1]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])

    # reorder K5
    mu_k5 <- rep(0, K[5])
    for(i in 1:K[5]) {
      mu_k5[i] <- mu_matrix[1, 1, 1, 1, i, 1] - mu_matrix[1, 1, 1, 1, ref[5], 1]
    }
    mu_k5_sort <- sort(mu_k5)
    mu <- c(mu, mu_k5_sort[mu_k5_sort != 0])
    # order of re-arranged cluster for omics 5
    k5 <- order(mu_k5)
    k5_order <- c(ref[5], k5[k5 != ref[5]])

    # reorder K6
    mu_k6 <- rep(0, K[6])
    for(i in 1:K[6]) {
      mu_k6[i] <- mu_matrix[1, 1, 1, 1, 1, i] - mu_matrix[1, 1, 1, 1, 1, ref[6]]
    }
    mu_k6_sort <- sort(mu_k6)
    mu <- c(mu, mu_k6_sort[mu_k6_sort != 0])
    # order of re-arranged cluster for omics 6
    k6 <- order(mu_k6)
    k6_order <- c(ref[6], k6[k6 != ref[6]])

    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order,
                    K5 = k5_order,
                    K6 = k6_order)
  }


  Delta$mu <- mu
  return(list(Delta = Delta,
              K_order = K_order))
}



reorder_Mu_Sigma <- function(Mu_Sigma, K_order) {
  for(i in 1:length(K_order)) {
    temp_Mu <- Mu_Sigma$Mu[[i]]
    temp_Sigma <- Mu_Sigma$Sigma[[i]]
    # reorder Mu
    Mu_Sigma$Mu[[i]] <- temp_Mu[, K_order[[i]]]
    Mu_Sigma$Sigma[[i]] <- temp_Sigma[, , K_order[[i]]]
  }
  return(Mu_Sigma)
}



reorder_Beta <- function(Beta, K_order) {
  for(i in 1:length(K_order)) {
    temp_Beta <- Beta[[i]]
    temp_Beta <- rbind(rep(0, ncol(temp_Beta)),
                       temp_Beta)
    temp_Beta_reorder <- temp_Beta[K_order[[i]], ]
    ref <- temp_Beta_reorder[1, ]
    for(j in 1:nrow(temp_Beta_reorder)) {
      temp_Beta_reorder[j, ] <- temp_Beta_reorder[j, ] - ref
    }
    Beta[[i]] <- temp_Beta_reorder[-1, ]
  }

  return(Beta)
}



reorder_z <- function(z, K_order) {
  if(length(K_order) == 2) {
    z <- z[K_order[[1]], K_order[[2]], ]
  }
  return(z)
}


###function to reorder all model parameters###
reorder_lucid <- function(model) {
  ref <- get_ref_cluster(Delta = model$res_Delta$Delta)
  r_Delta <- reorder_Delta(ref = ref,
                           Delta = model$res_Delta$Delta)
  r_Mu_Sigma <- reorder_Mu_Sigma(model$res_Mu_Sigma,
                                 K_order = r_Delta$K_order)
  r_Beta <- reorder_Beta(Beta = model$res_Beta$Beta,
                         K_order = r_Delta$K_order)
  model$res_Delta$Delta <- r_Delta
  model$res_Mu_Sigma$Mu <- r_Mu_Sigma$Mu
  model$res_Mu_Sigma$Sigma <- r_Mu_Sigma$Sigma
  model$res_Beta$Beta <- r_Beta
  model$z <- reorder_z(model$z, K_order = r_Delta$K_order)
  return(model)
}


# function to calculate BIC for LUCID in parallel
cal_bic_parallel <- function(object) {
  ##not having regularity yet, to be added
  s1 <- object$select$selectG
  s2 <- object$select$selectZ
  nG <- sum(s1)
  nZ <- sapply(s2,sum)
  K <- object$K
  #obtain number of parameters
  if(object$family == "gaussian"){
    nY <- length(object$res_Gamma$Gamma$mu) + length(object$res_Gamma$Gamma$sd)
  }
  if(object$family == "binomial"){
    #binary summary res_Gamma$Gamma$mu is unclear, use object$res_Gamma$fit$coefficients instead
    nY <- length(object$res_Gamma$fit$coefficients)
  }
  #initiate number of parameters
  npars = 0
  #compute number of parameters for G to X association
  for (i in 1:length(K)){
    npars_new = (nG + 1) * (K[i] - 1)
    npars = npars + npars_new
  }
  #compute number of parameters for X to Z association and add
  for (i in 1:length(K)){
    npars_new = (nZ[i] * K[i] + nZ[i] * nZ[i] * K[i])
    npars = npars + npars_new
  }
  #compute number of parameters for X to Y association and add
  npars = npars + nY

  BIC <- -2 * object$likelihood + npars * log(nrow(object$inclusion.p[[1]]))

  return(BIC)
}

# function to calculate BIC for LUCID in serial
cal_bic_serial <- function(object) {

  sum_all_sub = summary_lucid_simple(object)
  BIC = 0
  for (i in 1:length(sum_all_sub)){
    BIC_temp = sum_all_sub[[i]]$BIC
    BIC = BIC + BIC_temp
    }
  return(BIC)
}

cal_loglik_serial <- function(object) {

  sum_all_sub = summary_lucid_simple(object)
  loglik = 0
  for (i in 1:length(sum_all_sub)){
    loglik_temp = sum_all_sub[[i]]$loglik
    loglik = loglik + loglik_temp
  }
  return(loglik)
}


########Auxilary functions for printing the summary for LUCID in serial.############

print.auxi.serial.scen1 <- function(x, ...){
  if(inherits(x, "sumlucid_early")){
    K <- x$K
    beta <- as.data.frame(x$beta)
    dim1 <- ncol(beta) - 1
    z.mean <- as.data.frame(t(x$mu))
    cat("----------Summary of the First Component of the LUCID in Serial Model---------- \n \n")
    cat("----------Summary of the LUCID Early Integration Sub Model---------- \n \n")
    cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")


    cat("(1) Z: mean of omics data for each latent cluster \n")
    if(is.null(x$boot.se)){
      colnames(z.mean) <- paste0("mu_cluster", 1:K)
      print(z.mean)
    } else{
      print(x$boot.se$mu)
    }
    cat("\n")
    cat("(2) E: odds ratio of being assigned to each latent cluster for each exposure \n")
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
    cat("----------Summary of the First Component of the LUCID in Serial Model---------- \n \n")
    cat("----------Summary of the LUCID in Parallel Sub Model---------- \n \n")
    cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")

    cat("(1) Z: mean of omics data for each latent cluster of each layer \n")
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
    cat("(2) E: odds ratio of being assigned to each latent cluster for each exposure for each layer \n")
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


print.auxi.serial.scen2 <- function(x, ...){
  if(inherits(x, "sumlucid_early")){
    K <- x$K
    delta <- as.data.frame(x$beta)
    dim1 <- ncol(delta) - 1
    z.mean <- as.data.frame(t(x$mu))
    cat("----------Summary of the Middle Component of the LUCID in Serial Model---------- \n \n")
    cat("----------Summary of the LUCID Early Integration Sub Model---------- \n \n")
    cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")


    cat("(1) Z: mean of omics data for each latent cluster \n")
    if(is.null(x$boot.se)){
      colnames(z.mean) <- paste0("mu_cluster", 1:K)
      print(z.mean)
    } else{
      print(x$boot.se$mu)
    }
    cat("\n")
    cat("(2) E: odds ratio of being assigned to each latent cluster for each cluster from the last sub model \n")
    if(is.null(ncol(delta))) {
      cat("no cluster is selected given current penalty Rho_G, please use a smaller penalty")
    } else {
      dd <- as.matrix(as.data.frame(delta)[2:K, 2:ncol(delta)])
      g.or <- data.frame(delta = unlist(split(dd, row(dd))))
      rownames(g.or) <- paste0(colnames(delta)[-1], ".cluster", sapply(2:K, function(x) return(rep(x, dim1))))
      if(is.null(x$boot.se)){
        g.or$OR <- exp(g.or$delta)
        print(g.or)
      } else{
        print(x$boot.se$delta)
      }
    }
  }
  if(inherits(x, "sumlucid_parallel")){
    K <- x$K
    #list of deltas for each layer
    delta <- x$beta$Beta
    dim1 <- ncol(x$beta$Beta[[1]]) - 1
    #list of Z means for each layer, transposed
    z.mean <- x$mu
    cat("----------Summary of the Middle Component of the LUCID in Serial Model---------- \n \n")
    cat("----------Summary of the LUCID in Parallel Sub Model---------- \n \n")
    cat("K = ", K, ", log likelihood =", x$loglik, ", BIC = ", x$BIC, "\n \n")

    cat("(1) Z: mean of omics data for each latent cluster of each layer \n")
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
    cat("(2) E: odds ratio of being assigned to each latent cluster for each cluster from the last sub model for each layer \n")
    if(is.null(delta)) {
      cat("no cluster is selected given current penalty Rho_G, please use a smaller penalty")
    } else {
      for (i in 1:length(K)){
        cat("Layer ",i, "\n \n")
        dd <- as.matrix(as.data.frame(delta[[i]][, -1]))
        g.or <- data.frame(delta = unlist(split(dd, row(dd))))

        rownames(g.or) <- paste0(colnames(delta[[i]])[-1], ".cluster", sapply(2:K[i], function(x) return(rep(x, dim1))))


        if(is.null(x$boot.se)){
          g.or$OR <- exp(g.or$delta)
          print(g.or)
        } else{
          print(x$boot.se$delta)
        }
        cat("\n \n")
      }
    }
  }
}


print.auxi.serial.scen3 <- function(x, ...){
  if(inherits(x, "sumlucid_early")){
    K <- x$K
    delta <- as.data.frame(x$beta)
    dim1 <- ncol(delta) - 1
    z.mean <- as.data.frame(t(x$mu))
    cat("----------Summary of the Last Component of the LUCID in Serial Model---------- \n \n")
    cat("----------Summary of the LUCID Early Integration Sub Model---------- \n \n")
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
    cat("(3) E: odds ratio of being assigned to each latent cluster for each cluster from the last sub model \n")
    if(is.null(ncol(delta))) {
      cat("no cluster is selected given current penalty Rho_G, please use a smaller penalty")
    } else {
      dd <- as.matrix(as.data.frame(delta)[2:K, 2:ncol(delta)])
      g.or <- data.frame(delta = unlist(split(dd, row(dd))))
      rownames(g.or) <- paste0(colnames(delta)[-1], ".cluster", sapply(2:K, function(x) return(rep(x, dim1))))
      if(is.null(x$boot.se)){
        g.or$OR <- exp(g.or$delta)
        print(g.or)
      } else{
        print(x$boot.se$delta)
      }
    }
  }
  if(inherits(x, "sumlucid_parallel")){
    K <- x$K
    #list of deltas for each layer
    delta <- x$beta$Beta
    dim1 <- ncol(x$beta$Beta[[1]]) - 1
    #list of Z means for each layer, transposed
    z.mean <- x$mu
    cat("----------Summary of the Last Component of the LUCID in Serial Model---------- \n \n")
    cat("----------Summary of the LUCID in Parallel Sub Model---------- \n \n")
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
    cat("(3) E: odds ratio of being assigned to each latent cluster for each cluster from the last sub model for each layer \n")
    if(is.null(delta)) {
      cat("no cluster is selected given current penalty Rho_G, please use a smaller penalty")
    } else {
      for (i in 1:length(K)){
        cat("Layer ",i, "\n \n")
        dd <- as.matrix(as.data.frame(delta[[i]][, -1]))
        g.or <- data.frame(delta = unlist(split(dd, row(dd))))

        rownames(g.or) <- paste0(colnames(delta[[i]])[-1], ".cluster", sapply(2:K[i], function(x) return(rep(x, dim1))))


        if(is.null(x$boot.se)){
          g.or$OR <- exp(g.or$delta)
          print(g.or)
        } else{
          print(x$boot.se$delta)
        }
        cat("\n \n")
      }
    }
  }
}


summary_lucid_simple <- function(object, boot.se = NULL){
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

    #return a list of sumlucid object for each submodel
    class(summary.list) <- "sumlucid_serial"
    return(summary.list)
  }
}

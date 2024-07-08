# Calculate the log-likelihood of cluster assignment for each observation for g-com of early
Estep_early_g <- function(beta,
                          mu,
                          sigma,
                          gamma = NULL,
                          G,
                          Z,
                          Y = NULL,
                          family.list,
                          K,
                          N,
                          useY,
                          ind.na, ...) {
  # initialize vectors for storing likelihood
  pXgG <- pZgX <- pYgX <- matrix(rep(0, N * K), nrow = N)
  
  # log-likelihood for G -> X
  xb <- cbind(rep(1, N), G) %*% t(beta)
  xb_lse <- apply(xb, 1, lse)
  pXgG <- xb - xb_lse
  
  
  
  vec <- pXgG 
  return (vec)
}


# Calculate the log-likelihood of cluster assignment for each observation for g-com of parallel
Estep_g <- function(G, Z, Y, Beta, Mu, Sigma, Delta, family, useY, na_pattern) {
  N <- nrow(G)
  K <- Delta$K
  res <- array(rep(0, prod(K) * N),
               dim = c(K, N))
  
  
  # E step for 2 omics data
  if(length(K) == 2) {
    f1 <- lapply(1:2, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })

    
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        res[i, j, ] <- f1[[1]][, i] + f1[[2]][, j] 
      }
    }
  }
  
  # E step for 3 omics data
  if(length(K) == 3) {
    f1 <- lapply(1:3, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          res[i, j, k, ] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] 
        }
      }
    }
  }
  
  
  # E step for 4 omics layers
  if(length(K) == 4) {
    f1 <- lapply(1:4, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          for(l in 1:K[4]) {
            res[i, j, k, l, ] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] + f1[[4]][, l] 
          }
        }
      }
    }
  }
  
  
  
  # E step for 5 omics layers
  if(length(K) == 5) {
    f1 <- lapply(1:5, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })

    
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          for(l in 1:K[4]) {
            for(m in 1:K[5]) {
              res[i, j, k, l, m,] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] + f1[[4]][, l] + f1[[5]][, m] 
            }
          }
        }
      }
    }
  }
  
  
  return(res)
}





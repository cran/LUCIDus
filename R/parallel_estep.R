# parameters
# 1 - number of latent cluster: K
# 2 - G -> X coef: Beta
# 3 - X -> Z coef: Mu, Sigma
# 4 - X -> Y coef: Delta



#Calculate the log prior inclusion probability for omics layer X_j given the exposures G

f_GtoX <- function(G, Beta_matrix) {
  N <- nrow(G)
  Beta_matrix <- rbind(rep(0, ncol(Beta_matrix)),
                       Beta_matrix)
  xb <- cbind(rep(1, N), G) %*% t(Beta_matrix)
  xb_LSE <- apply(xb, 1, LogSumExp)
  return(xb - xb_LSE)
}


# calculate the log likelihood for omics layer Z_j given the latent cluster X_j

f_XtoZ <- function(Z, Mu_matrix, Sigma_matrix) {
  N <- nrow(Z)
  M <- ncol(Z)
  K <- nrow(Mu_matrix)  # Number of clusters
  
  # Validate dimensions
  if (ncol(Mu_matrix) != M) {
    stop(sprintf("Dimension mismatch: Z has %d features but Mu_matrix has %d features", 
                M, ncol(Mu_matrix)))
  }
  if (!identical(dim(Sigma_matrix)[1:2], c(M, M))) {
    stop(sprintf("Invalid Sigma dimensions: Expected first two dimensions to be c(%d,%d), got c(%s)", 
                M, M, paste(dim(Sigma_matrix)[1:2], collapse=",")))
  }
  if (dim(Sigma_matrix)[3] != K) {
    stop(sprintf("Invalid Sigma dimensions: Expected %d clusters, got %d", 
                K, dim(Sigma_matrix)[3]))
  }
  
  XtoZ <- matrix(rep(0, N * K), nrow = N)
  for (i in 1:K) {
    mean_i <- Mu_matrix[i, ]
    sigma_i <- Sigma_matrix[, , i]
    
    # Additional validation
    if (length(mean_i) != M) {
      stop(sprintf("Mean vector for cluster %d has wrong length: expected %d, got %d", 
                  i, M, length(mean_i)))
    }
    
    # Ensure sigma is symmetric and positive definite
    sigma_i <- (sigma_i + t(sigma_i))/2  # Ensure symmetry
    eigen_values <- eigen(sigma_i, symmetric = TRUE)$values
    if (any(eigen_values <= 0)) {
      # Add small constant to diagonal if not positive definite
      sigma_i <- sigma_i + diag(1e-6, M)
    }
    
    tmp_ll <- try({
      mclust::dmvnorm(data = Z,
                      mean = mean_i,
                      sigma = sigma_i,
                      log = TRUE)
    }, silent = TRUE)
    
    # Check for errors
    if (inherits(tmp_ll, "try-error")) {
      stop(sprintf("Error in cluster %d: %s", i, attr(tmp_ll, "condition")$message))
    }
    XtoZ[, i] <- tmp_ll
  }
  return(XtoZ)
}


# Calculate the log likelihood of outcome Y given all latent variables X

f_XtoY <- function(Y, Delta, family) {
  family <- to_parallel_family(family)

  if(!is.matrix(Y)) {
    stop("Y should be a matrix")
  }
  N <- nrow(Y)
  K <- Delta$K
  XtoY <- array(rep(0, prod(K) * N), dim = c(K, N))

  # if the outcome is continuous
  if(family == "gaussian") {

    # if 1 omics layers
    if(length(K) == 1) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for(i in 1:K[1]) {
        XtoY[i] <- dnorm(Y, mean = mu[i], sd = Delta$sd, log = TRUE)
        }
      }



    # if 2 omics layers
    if(length(K) == 2) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          XtoY[i, j, ] <- dnorm(Y, mean = mu[i, j], sd = Delta$sd, log = TRUE)
        }
      }
    }

    # if 3 omics layers
    if(length(K) == 3) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for (i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            XtoY[i, j, k, ] <- dnorm(Y, mean = mu[i, j, k], sd = Delta$sd, log = TRUE)
          }
        }
      }
    }

    # if 4 omics layers
    if(length(K) == 4) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for (i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            for(l in 1:K[4]) {
              XtoY[i, j, k, l, ] <- dnorm(Y, mean = mu[i, j, k, l], sd = Delta$sd, log = TRUE)
            }
          }
        }
      }
    }

    # if 5 omics layers
    if(length(K) == 5) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for (i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            for(l in 1:K[4]) {
              for(m in 1:K[5]) {
                XtoY[i, j, k, l, m, ] <- dnorm(Y, mean = mu[i, j, k, l, m], sd = Delta$sd, log = TRUE)
              }
            }
          }
        }
      }
    }


  }

  # if the outcome is binary
  if(family == "binomial") {
    p <- Delta$mu

    # if 1 omics layers
    if(length(K) == 1) {
      for(i in 1:K[1]) {
        XtoY[i] <- dbinom(Y, size = 1, prob = p[i], log = TRUE)
        }
      }


    # if 2 omics layers
    if(length(K) == 2) {
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          XtoY[i, j, ] <- dbinom(Y, size = 1, prob = p[i, j], log = TRUE)
        }
      }
    }

    # if 3 omics layers
    if(length(K) == 3) {
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            XtoY[i, j, k, ] <- dbinom(Y, size = 1, prob = p[i, j, k], log = TRUE)
          }
        }
      }
    }

    # if 4 omics layers
    if(length(K) == 4) {
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            for(l in 1:K[4]) {
              XtoY[i, j, k, l, ] <- dbinom(Y, size = 1, prob = p[i, j, k, l], log = TRUE)
            }
          }
        }
      }
    }

    # if 5 omics layers
    if(length(K) == 5) {
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            for(l in 1:K[4]) {
              for(m in 1:K[5]) {
                XtoY[i, j, k, l, m, ] <- dbinom(Y, size = 1, prob = p[i, j, k, l, m], log = TRUE)
              }
            }
          }
        }
      }
    }

  }
  return(XtoY)
}


# Estep
Estep <- function(G, Z, Y, Beta, Mu, Sigma, Delta, family, useY, na_pattern) {
  family <- to_parallel_family(family)
  N <- nrow(G)
  K <- Delta$K
  res <- array(rep(0, prod(K) * N),
               dim = c(K, N))
  
  # Build layer-wise Z|X log-likelihood matrices on the full N rows.
  # Rows with all-missing omics values contribute 0 (no Z information).
  f2_all <- lapply(seq_along(K), function(i) {
    layer_res <- matrix(0, nrow = N, ncol = K[i])
    idx_obs <- na_pattern[[i]]$indicator_na != 3
    if (any(idx_obs)) {
      layer_res[idx_obs, ] <- f_XtoZ(
        Z = Z[[i]][idx_obs, , drop = FALSE],
        Mu_matrix = Mu[[i]],
        Sigma_matrix = Sigma[[i]]
      )
    }
    layer_res
  })


  # E step for 2 omics data
  if(length(K) == 2) {
    f1 <- lapply(1:2, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    if(useY) {
      f3 <- f_XtoY(Y = Y, Delta = Delta, family = family)
    }

    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        res[i, j, ] <- f1[[1]][, i] + f1[[2]][, j] + f2_all[[1]][, i] + f2_all[[2]][, j]
        if(useY) {
          res[i, j, ] <- res[i, j, ] + f3[i, j, ]
        }
      }
    }
  }

  # E step for 3 omics data
  if(length(K) == 3) {
    f1 <- lapply(1:3, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    if(useY) {
      f3 <- f_XtoY(Y = Y, Delta = Delta, family = family)
    }

    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          res[i, j, k, ] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] + f2_all[[1]][, i] + f2_all[[2]][, j] + f2_all[[3]][, k]
          if(useY) {
            res[i, j, k, ] <- res[i, j, k, ] + f3[i, j, k, ]
          }
        }
      }
    }
  }


  # E step for 4 omics layers
  if(length(K) == 4) {
    f1 <- lapply(1:4, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    if(useY) {
      f3 <- f_XtoY(Y = Y, Delta = Delta, family = family)
    }

    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          for(l in 1:K[4]) {
            res[i, j, k, l, ] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] + f1[[4]][, l] + f2_all[[1]][, i] + f2_all[[2]][, j] + f2_all[[3]][, k] + f2_all[[4]][, l]
            if(useY) {
              res[i, j, k, l, ] <- res[i, j, k, l, ] + f3[i, j, k, l, ]
            }
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
    if(useY) {
      f3 <- f_XtoY(Y = Y, Delta = Delta, family = family)
    }

    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          for(l in 1:K[4]) {
            for(m in 1:K[5]) {
              res[i, j, k, l, m,] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] + f1[[4]][, l] + f1[[5]][, m] + f2_all[[1]][, i] + f2_all[[2]][, j] + f2_all[[3]][, k] + f2_all[[4]][, l] + f2_all[[5]][, m]
              if(useY) {
                res[i, j, k, l, m,] <- res[i, j, k, l, m, ] + f3[i, j, k, l, m, ]
              }
            }
          }
        }
      }
    }
  }


  return(res)
}



Estep_to_r <- function(Estep_array, K, N) {

  # E step for 1 omics layers
  if(length(K) == 1) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, i] <- array(r, dim = K)
    }
  }

  # E step for 2 omics layers
  if(length(K) == 2) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, , i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, , i] <- array(r, dim = K)
    }
  }

  # E step for 3 omics layers
  if(length(K) == 3) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, , , i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, , , i] <- array(r, dim = K)
    }
  }


  # E step for 4 omics layers
  if(length(K) == 4) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, , , , i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, , , , i] <- array(r, dim = K)
    }
  }

  # E step for 5 omics layers
  if(length(K) == 5) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, , , , , i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, , , , , i] <- array(r, dim = K)
    }
  }


  return(Estep_array)
}

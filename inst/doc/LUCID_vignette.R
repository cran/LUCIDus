## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(LUCIDus)

## ----out.width="50%", echo=FALSE----------------------------------------------
knitr::include_graphics("DAG.png")

## ---- eval=FALSE--------------------------------------------------------------
#  library(LUCIDus)
#  # use simulated data
#  G <- sim_data$G
#  Z <- sim_data$Z
#  Y_normal <- sim_data$Y_normal
#  Y_binary <- sim_data$Y_binary
#  cov <- sim_data$Covariate
#  
#  # fit LUCID model with continuous outcome
#  fit1 <- est.lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, seed = 1008)
#  
#  
#  # fit LUCID model with binary outcome
#  fit2 <- est.lucid(G = G, Z = Z, Y = Y_binary, family = "binary", K = 2, seed = 1008)
#  
#  # fit LUCID model with covariates
#  fit3 <- est.lucid(G = G, Z = Z, Y = Y_binary, CoY = cov, family = "binary", K = 2, seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # fit LUCID model without useing information from outcome
#  fit4 <- est.lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, useY = FALSE, seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # fit LUCID model with automatic selection on optimal covariance models
#  fit5 <- est.lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, modelName = NULL, seed = 1008)
#  # check the optimal model
#  fit5$modelName
#  
#  # fit LUCID model with a specified covariance model
#  fit6 <- est.lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, modelName = "EII", seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # initialize EM algorithm by mclust
#  fit7 <- est.lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, init_par = "mclust" , seed = 1008)
#  
#  # initialize EM algorithm via randomization
#  fit8 <- est.lucid(G = G, Z = Z, Y = Y_normal, family = "normal", K = 2, init_par = "random" , seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # summarize lucid results
#  summary_lucid(fit1)

## ---- eval=FALSE--------------------------------------------------------------
#  # visualze lucid model via a Snakey diagram
#  plot_lucid(fit1)

## ----out.width="50%", include=FALSE-------------------------------------------
knitr::include_graphics("sankey.png")

## ---- eval=FALSE--------------------------------------------------------------
#  # change node color
#  plot_lucid(fit1, G_color = "yellow")
#  plot_lucid(fit1, Z_color = "red")
#  
#  # change link color
#  plot_lucid(fit1, pos_link_color = "red", neg_link_color = "green")

## ---- eval=FALSE--------------------------------------------------------------
#  # fit LUCID model with block-wise missing pattern in omics data
#  Z_miss_1 <- Z
#  Z_miss_1[sample(1:nrow(Z), 0.3 * nrow(Z)), ] <- NA
#  fit9 <- est.lucid(G = G, Z = Z_miss_1, Y = Y_normal, family = "normal", K = 2)
#  
#  # fit LUCID model with sporadic missing pattern in omics data
#  Z_miss_2 <- Z
#  index <- arrayInd(sample(length(Z_miss_2), 0.3 * length(Z_miss_2)), dim(Z_miss_2))
#  Z_miss_2[index] <- NA
#  fit10 <- est.lucid(G = G, Z = Z_miss_2, Y = Y_normal, family = "normal", K = 2, seed = 1008)

## ---- eval=FALSE--------------------------------------------------------------
#  # imputed omics dataset
#  fit10$Z

## ---- eval=FALSE--------------------------------------------------------------
#  # use LUCID model to conduct integrated variable selection
#  # select exposure
#  fit6 <- est.lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal",
#  K = 2, seed = 1008, Rho_G = 0.1)
#  # select omics data
#  fit7 <- est.lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal",
#  K = 2, seed = 1008, Rho_Z_Mu = 90, Rho_Z_Cov = 0.1, init_par = "random")

## ---- eval=FALSE--------------------------------------------------------------
#  # tune lucid over a grid of K (note this function may take time to run)
#  tune_lucid <- lucid(G = G, Z = Z, Y = Y_normal, K =2:5)

## ---- eval=FALSE--------------------------------------------------------------
#  > tune_lucid$tune_list
#    K Rho_G Rho_Z_Mu Rho_Z_Cov      BIC
#  1 2     0        0         0 46917.94
#  2 3     0        0         0 47716.94
#  3 4     0        0         0 48520.86
#  4 5     0        0         0 49276.93

## ---- eval=FALSE--------------------------------------------------------------
#  tune_lucid$best_model

## ---- eval=FALSE--------------------------------------------------------------
#  # conduct bootstrap resampling
#  boot1 <- boot.lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100)
#  
#  # use 90% CI
#  boot2 <- boot.lucid(G = G, Z = Z, Y = Y_normal, model = fit1, R = 100, conf = 0.9)

## ---- eval=FALSE--------------------------------------------------------------
#  # check distribution for bootstrap replicates of the variable of interest
#  plot(boot1$bootstrap, 1)


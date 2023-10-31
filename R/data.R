#' @title A simulated dataset for LUCID 
#'
#' @description This is an example dataset to illustrate LUCID model. It is simulated
#' by assuming there are 2 latent clusters in the data. We assume the exposures 
#' are associated with latent cluster which ultimately affects the PFAS concentration
#' and liver injury in children. The latent clusters are also characterized by 
#' differential levels of metabolites.
#' 
#' @format A list with 5 matrices corresponding to exposures (G), omics data (Z),
#' a continuous outcome, a binary outcome and 2 covariates (can be used either 
#' as CoX or CoY). Each matrice contains 2000 observations.
#' \describe{
#'     \item{G}{10 exposures}
#'     \item{Z}{10 metabolites}
#'     \item{Y_normal}{Outcome, PFAS concentration in children}
#'     \item{Y_bninary}{Bianry outcome, liver injury status}
#'     \item{Covariates}{2 continous covariates, can be treated as either CoX or 
#'                       CoY}
#'     \item{X}{Latent clusters}
#' }
"sim_data"


#' @title A simulated HELIX dataset for LUCID  
#'
#' @description The Human Early-Life Exposome (HELIX) project is multi-center 
#' research project that aims to characterize early-life environmental exposures
#' and associate these with omics biomarkers and child health outcomes (Vrijheid, 2014. doi: 10.1289/ehp.1307204).
#' We used a subset of HELIX data from Exposome Data Challenge 2021 (hold by
#' ISGlobal) as an example to illustrate LUCID model.
#' 
#' @format A list with 4 matrices corresponding to exposures (G), omics data (Z),
#' outcome (Y) and covariates (CoY), a total of 420 observations
#' \describe{
#'     \item{exposure}{1 exposures measuring the maternal exposure to utero mercury.}
#'     \item{omics}{10 methylomics, 10 transcriptomics, 10 miRNA}
#'     \item{outcome}{A continuous outcome as an indicator of metabolic-dysfunciton-associated fatty liver
#'     disease (MAFLD)}
#'     \item{covariate}{3 covariates including fish_preg_ter, child sex, maternal
#'     age}
#' }
"simulated_HELIX_data"
#' @title simulated dataset 1
#'
#' @description A simulated dataset for integrated clustering with normal outcome. The data is simulated under cluster number K = 2.
#' @format A matrix of 22 columns, which are
#' \describe{
#'   \item{G1 - G10}{Genetic features, G1 to G5 are causal genes contributed to clustering, with OR = 2; G6 to G10 are null genes that is not related to clustering}
#'   \item{Z1 - Z10}{Biomarkers, Z1 to Z5 are causal biomarkers with delta Z = 4 between 2 clusters, Z6 to Z10 are noises with delta Z = 0. All biomarkers are assumed to be independent with each other}
#'   \item{Y}{Outcome of interest, which follows 2 normal distribution with N(-1, 1) and N(1, 1)}
#'   \item{X}{Latent cluster assignment for each observation}
#' }
"sim1"

#' @title simulated dataset 2
#'
#' @description A simulated dataset for integrated clustering with binary outcome. The data is simulated under cluster number K = 2.
#' @format A matrix of 22 columns, which are
#' \describe{
#'   \item{G1 - G10}{Genetic features, G1 to G5 are causal genes contributed to clustering, with OR = 2; G6 to G10 are null genes that is not related to clustering}
#'   \item{Z1 - Z10}{Biomarkers, Z1 to Z5 are causal biomarkers with delta Z = 4 between 2 clusters, Z6 to Z10 are noises with delta Z = 0. All biomarkers are assumed to be independent with each other}
#'   \item{Y}{Outcome of interest, the odds ratio of the cluster is 2}
#'   \item{X}{Latent cluster assignment for each observation}
#' }
"sim2"

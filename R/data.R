#' @title A simulated dataset for LUCID
#'
#' @description An example dataset used to illustrate the LUCID model, simulated
#' under two latent clusters. The exposures are associated with the latent
#' cluster, which in turn affects PFAS concentration and liver injury in
#' children; the clusters are also characterised by differential metabolite
#' levels. Because the generating cluster membership is retained in \code{X}, the
#' data can be used to check recovered clusters against the truth, up to the
#' arbitrary labelling of clusters.
#'
#' @format A list of 6 elements, each with 2000 observations:
#' \describe{
#'     \item{G}{A 2000 by 10 matrix of exposures.}
#'     \item{Z}{A 2000 by 10 matrix of metabolites.}
#'     \item{Y_normal}{A 2000 by 1 matrix; continuous outcome, PFAS
#'                     concentration in children.}
#'     \item{Y_binary}{A 2000 by 1 matrix; binary outcome, liver injury status,
#'                     coded 0 and 1.}
#'     \item{Covariate}{A 2000 by 2 matrix of continuous covariates, usable as
#'                      either \code{CoG} or \code{CoY}.}
#'     \item{X}{An integer vector of length 2000 giving the latent cluster each
#'              observation was generated from.}
#' }
"sim_data"


#' @title A simulated HELIX dataset for LUCID
#'
#' @description The Human Early-Life Exposome (HELIX) project is a multi-center
#' research project that aims to characterize early-life environmental exposures
#' and associate these with omics biomarkers and child health outcomes (Vrijheid,
#' 2014. doi: 10.1289/ehp.1307204). This is a simulated subset of the HELIX data
#' released for the Exposome Data Challenge 2021 (held by ISGlobal), used to
#' illustrate the LUCID model on three omics layers.
#'
#' The three omics layers share the same 420 subjects in the same row order, so
#' they can be passed together as the \code{Z} list of a parallel or serial fit,
#' or used one at a time for early integration.
#'
#' @format A list of 4 elements, each with 420 observations:
#' \describe{
#'     \item{phenotype}{A 420 by 6 data frame holding the exposure, the outcome
#'       and the covariates: \code{id}; \code{hs_hg_m_scaled}, the scaled
#'       maternal exposure to in-utero mercury, used as \code{G};
#'       \code{ck18_scaled}, a scaled continuous indicator of
#'       metabolic-dysfunction-associated fatty liver disease (MAFLD), used as
#'       \code{Y}; and the covariates \code{hs_child_age_yrs_None},
#'       \code{e3_sex_None} and \code{h_fish_preg_Ter}.}
#'     \item{methylome}{A 420 by 10 matrix of methylomics features.}
#'     \item{transcriptome}{A 420 by 10 matrix of transcriptomics features.}
#'     \item{miRNA}{A 420 by 10 matrix of miRNA features.}
#' }
"simulated_HELIX_data"

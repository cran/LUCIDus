% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missing.R
\name{check_na}
\alias{check_na}
\title{Check missing patterns in one layer of omics data Z}
\usage{
check_na(Z)
}
\arguments{
\item{Z}{A data matrix representing one layer of omics data}
}
\value{
\enumerate{
\item index:indeces for missing values in the omics data
\item indicator_na: missing pattern for each observation
\item impute_flag: - flag to initialize imputation. Only happens when sporadic missing
pattern is observed
}
}
\description{
Check missing patterns in one layer of omics data Z
}

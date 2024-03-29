% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missing.R
\name{fill_data}
\alias{fill_data}
\title{Impute missing data by optimizing the likelihood function}
\usage{
fill_data(obs, mu, sigma, p, index, lucid_model)
}
\arguments{
\item{obs}{a vector of length M}

\item{mu}{a matrix of size M x K}

\item{sigma}{a matrix of size M x M x K}

\item{p}{a vector of length K}

\item{index}{a vector of length M, indicating whether a value is missing
or not in the raw data}

\item{lucid_model}{Specifying LUCID model, "early" for early integration, "parallel" for lucid in parallel,
"serial" for lucid in serial}
}
\value{
an observation with updated imputed value
}
\description{
Impute missing data by optimizing the likelihood function
}

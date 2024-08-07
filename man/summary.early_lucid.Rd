% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.early_lucid}
\alias{summary.early_lucid}
\title{Summarize results of the early LUCID model}
\usage{
\method{summary}{early_lucid}(object, ...)
}
\arguments{
\item{object}{A LUCID model fitted by \code{\link{estimate_lucid}}}

\item{...}{Additional argument of \code{boot.se}, which is an object returned by \code{\link{boot_lucid}}}
}
\value{
A list, containing the extracted key parameters from the LUCID model that can be used to print the summary table
}
\description{
Summarize results of the early LUCID model
}
\examples{
\donttest{
# use simulated data
G <- sim_data$G
Z <- sim_data$Z
Y_normal <- sim_data$Y_normal

# fit lucid model
fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", family = "normal", K = 2,
seed = 1008)

# conduct bootstrap resampling
boot1 <- boot_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", model = fit1, R = 100)

# summarize lucid model
summary(fit1)

# summarize lucid model with bootstrap CIs
summary(fit1, boot.se = boot1)
}
}

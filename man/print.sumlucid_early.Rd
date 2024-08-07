% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{print.sumlucid_early}
\alias{print.sumlucid_early}
\title{Print the output of LUCID in a nicer table}
\usage{
\method{print}{sumlucid_early}(x, ...)
}
\arguments{
\item{x}{An object returned by \code{summary}}

\item{...}{Other parameters to be passed to \code{print.sumlucid_serial}}
}
\value{
A nice table/several nice tables of the summary of the LUCID model
}
\description{
Print the output of LUCID in a nicer table
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

# print the summary of the lucid model in a table
temp <- summary(fit1)
print(temp)
}
}

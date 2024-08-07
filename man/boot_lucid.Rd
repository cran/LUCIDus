% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_lucid.R
\name{boot_lucid}
\alias{boot_lucid}
\title{Inference of LUCID model based on bootstrap resampling}
\usage{
boot_lucid(
  G,
  Z,
  Y,
  lucid_model = c("early", "parallel", "serial"),
  CoG = NULL,
  CoY = NULL,
  model,
  conf = 0.95,
  R = 100,
  verbose = FALSE
)
}
\arguments{
\item{G}{Exposures, a numeric vector, matrix, or data frame. Categorical variable
should be transformed into dummy variables. If a matrix or data frame, rows
represent observations and columns correspond to variables.}

\item{Z}{Omics data for LUCID early integration, a numeric matrix or data frame. Rows correspond to observations
and columns correspond to variables.}

\item{Y}{Outcome, a numeric vector. Categorical variable is not allowed. Binary
outcome should be coded as 0 and 1.}

\item{lucid_model}{Specifying LUCID model, "early" for early integration, "parallel" for lucid in parallel,
"serial" for LUCID in serial.Now only work for LUCID early.
If "parallel" or "serial", the function will do nothing.}

\item{CoG}{Optional, covariates to be adjusted for estimating the latent cluster.
A numeric vector, matrix or data frame. Categorical variable should be transformed
into dummy variables.}

\item{CoY}{Optional, covariates to be adjusted for estimating the association
between latent cluster and the outcome. A numeric vector, matrix or data frame.
Categorical variable should be transformed into dummy variables.}

\item{model}{A LUCID model fitted by \code{estimate_lucid}.}

\item{conf}{A numeric scalar between 0 and 1 to specify confidence level(s)
of the required interval(s).}

\item{R}{An integer to specify number of bootstrap replicates for LUCID model.
If feasible, it is recommended to set R >= 1000.}

\item{verbose}{A flag indicates whether detailed information
is printed in console. Default is FALSE.}
}
\value{
A list, containing the following components:
\item{beta}{effect estimate for each exposure}
\item{mu}{cluster-specific mean for each omics feature}
\item{gamma}{effect estiamte for the association btween latent cluster and
outcome}
\item{bootstrap}{The \code{boot} object returned by \code{boot:boot}}
}
\description{
Generate \code{R} bootstrap replicates of LUCID parameters and
derive confidence interval (CI) base on bootstrap. Bootstrap replicates are
generated based on nonparameteric resampling, implemented by \code{ordinary}
method of \code{boot::boot} function. Now only achieved for LUCID early integration.
}
\examples{
\donttest{
# use simulated data
G <- sim_data$G
Z <- sim_data$Z
Y_normal <- sim_data$Y_normal

# fit lucid model
fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", 
family = "normal", K = 2,
seed = 1008)

# conduct bootstrap resampling
boot1 <- boot_lucid(G = G, Z = Z, Y = Y_normal, 
lucid_model = "early",model = fit1, R = 100)

# use 90\% CI
boot2 <- boot_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", 
model = fit1, R = 100, conf = 0.9)
}
}

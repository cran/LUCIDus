% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_lucid.R
\name{plot}
\alias{plot}
\title{Visualize LUCID model through a Sankey diagram}
\usage{
plot(x, ...)
}
\arguments{
\item{x}{A LUCID model fitted by \code{\link{estimate_lucid}}}

\item{...}{Additional arguments to specify colors and fontsize}
}
\value{
A DAG graph created by \code{\link[networkD3]{sankeyNetwork}}
}
\description{
In the Sankey diagram, each node either represents a variable (exposure,
omics or outcome) or a latent cluster. Each line represents an association. The
color of the node represents variable type, either exposure, omics or outcome.
The width of the line represents the effect size of a certain association; the
color of the line represents the direction of a certain association. Only work for LUCID early for now.
}
\examples{
# prepare data
G <- sim_data$G
Z <- sim_data$Z
Y_normal <- sim_data$Y_normal
Y_binary <- sim_data$Y_binary
cov <- sim_data$Covariate

# plot lucid model
fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", 
CoY = NULL, family = "normal", K = 2, seed = 1008)
plot(fit1)

# change node color
plot(fit1, G_color = "yellow")
plot(fit1, Z_color = "red")

# change link color
plot(fit1, pos_link_color = "red", neg_link_color = "green")
}

#' @title Visualize LUCID model through a Sankey diagram
#' @description In the Sankey diagram, each node either represents a variable (exposure,
#' omics or outcome) or a latent cluster. Each line represents an association. The
#' color of the node represents variable type, either exposure, omics or outcome.
#' The width of the line represents the effect size of a certain association; the
#' color of the line represents the direction of a certain association. Only work for LUCID early for now.
#'
#' @param x A LUCID model fitted by \code{\link{estimate_lucid}}
#' @param ... Additional arguments to specify colors and fontsize
#'
#' @return A DAG graph created by \code{\link[networkD3]{sankeyNetwork}}
#'
#' @import networkD3
#' @importFrom jsonlite toJSON
#'
#' @export
#'
#' @examples
#' # prepare data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#' Y_binary <- sim_data$Y_binary
#' cov <- sim_data$Covariate
#'
#' # plot lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early", 
#' CoY = NULL, family = "normal", K = 2, seed = 1008)
#' plot(fit1)
#'
#' # change node color
#' plot(fit1, G_color = "yellow")
#' plot(fit1, Z_color = "red")
#'
#' # change link color
#' plot(fit1, pos_link_color = "red", neg_link_color = "green")

# Define the generic plot function
#' @export
plot <- function(x, ...) {
  UseMethod("plot")
}

# Utility function to provide default values
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# Define plot.early_lucid function
#' @export
plot.early_lucid <- function(x, ...) {
  args <- list(...)
  G_color <- args$G_color %||% "dimgray"
  X_color <- args$X_color %||% "#eb8c30"
  Z_color <- args$Z_color %||% "#2fa4da"
  Y_color <- args$Y_color %||% "#afa58e"
  pos_link_color <- args$pos_link_color %||% "#67928b"
  neg_link_color <- args$neg_link_color %||% "#d1e5eb"
  fontsize <- args$fontsize %||% 7
  
  K <- x$K
  var.names <- x$var.names
  pars <- x$pars
  dimG <- length(var.names$Gnames)
  dimZ <- length(var.names$Znames)
  valueGtoX <- as.vector(t(x$res_Beta[, -1]))
  valueXtoZ <- as.vector(t(x$res_Mu))
  valueXtoY <- as.vector(x$res_Gamma$beta)[1:K]
  valueXtoY[1] <- 0
  GtoX <- data.frame(source = rep(x$var.names$Gnames, K),
                     target = paste0("Latent Cluster", as.vector(sapply(1:K, function(x) rep(x, dimG)))),
                     value = abs(valueGtoX),
                     group = as.factor(valueGtoX > 0))
  XtoZ <- data.frame(source = paste0("Latent Cluster", as.vector(sapply(1:K, function(x) rep(x, dimZ)))),
                     target = rep(var.names$Znames, K),
                     value = abs(valueXtoZ),
                     group = as.factor(valueXtoZ > 0))
  XtoY <- data.frame(source = paste0("Latent Cluster", 1:K),
                     target = rep(var.names$Ynames, K),
                     value = abs(valueXtoY),
                     group = as.factor(valueXtoY > 0))
  links <- rbind(GtoX, XtoZ, XtoY)
  nodes <- data.frame(name = unique(c(as.character(links$source), as.character(links$target))),
                      group = as.factor(c(rep("exposure", dimG),
                                          rep("lc", K),
                                          rep("biomarker", dimZ), "outcome")))
  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1
  color_scale <- data.frame(domain = c("exposure", "lc", "biomarker", "outcome", "TRUE", "FALSE"),
                            range = c(G_color, X_color, Z_color, Y_color, pos_link_color, neg_link_color))
  
  p <- sankeyNetwork(Links = links,
                     Nodes = nodes,
                     Source = "IDsource",
                     Target = "IDtarget",
                     Value = "value",
                     NodeID = "name",
                     colourScale = JS(
                       sprintf(
                         'd3.scaleOrdinal()
                               .domain(%s)
                               .range(%s)',
                         jsonlite::toJSON(color_scale$domain),
                         jsonlite::toJSON(color_scale$range)
                       )),
                     LinkGroup = "group",
                     NodeGroup = "group",
                     sinksRight = FALSE,
                     fontSize = fontsize)
  p
}

# Define plot.lucid_serial function
#' @export
plot.lucid_serial <- function(x, ...) {
  args <- list(...)
  G_color <- args$G_color %||% "dimgray"
  X_color <- args$X_color %||% "#eb8c30"
  Z_color <- args$Z_color %||% "#2fa4da"
  Y_color <- args$Y_color %||% "#afa58e"
  pos_link_color <- args$pos_link_color %||% "#67928b"
  neg_link_color <- args$neg_link_color %||% "#d1e5eb"
  fontsize <- args$fontsize %||% 7
  
  stop("The plotting function of LUCID in Parallel and Serial is still under development")
}

# Define plot.lucid_parallel function
#' @export
plot.lucid_parallel <- function(x, ...) {
  args <- list(...)
  G_color <- args$G_color %||% "dimgray"
  X_color <- args$X_color %||% "#eb8c30"
  Z_color <- args$Z_color %||% "#2fa4da"
  Y_color <- args$Y_color %||% "#afa58e"
  pos_link_color <- args$pos_link_color %||% "#67928b"
  neg_link_color <- args$neg_link_color %||% "#d1e5eb"
  fontsize <- args$fontsize %||% 7
  
  stop("The plotting function of LUCID in Parallel and Serial is still under development")
}

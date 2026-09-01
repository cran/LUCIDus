#' \code{NULL}-coalescing operator
#'
#' @param a A value, possibly \code{NULL}.
#' @param b Fallback used when \code{a} is \code{NULL}.
#' @return \code{a} if it is not \code{NULL}, otherwise \code{b}.
#' @noRd
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

#' Visualize an early-integration LUCID model through a Sankey diagram
#' @description
#' Draws the fitted model as a Sankey diagram: exposures flow into the latent
#' clusters, and the clusters flow on into the omics features and the outcome.
#' Each node is either a variable (exposure, omics or outcome) or a latent
#' cluster, and its colour indicates which. Each link is an estimated
#' association: its width is the magnitude of the effect and its colour the
#' sign, so the diagram shows at a glance which exposures drive which cluster
#' and how that cluster differs in the omics layer.
#'
#' Only exposures and omics features retained by the model are drawn, so a
#' penalized fit yields a correspondingly sparser diagram.
#'
#' @section Model types:
#' Implemented for early integration only. A \code{lucid_parallel} or
#' \code{lucid_serial} fit has a registered method
#' (\code{\link{plot.lucid_parallel}}/\code{\link{plot.lucid_serial}}), but
#' it raises an error: the parallel and serial diagrams are still under
#' development.
#'
#' @param x A LUCID model fitted by \code{\link{estimate_lucid}} or
#'   \code{\link{lucid}}, of class \code{early_lucid}.
#' @param ... Appearance options, all optional:
#'   \describe{
#'     \item{G_color}{Colour of the exposure nodes (default \code{"dimgray"}).}
#'     \item{X_color}{Colour of the latent-cluster nodes (default
#'       \code{"#eb8c30"}).}
#'     \item{Z_color}{Colour of the omics nodes (default \code{"#2fa4da"}).}
#'     \item{Y_color}{Colour of the outcome node (default \code{"#afa58e"}).}
#'     \item{pos_link_color}{Colour of links with a positive effect (default
#'       \code{"#67928b"}).}
#'     \item{neg_link_color}{Colour of links with a negative effect (default
#'       \code{"#d1e5eb"}).}
#'     \item{fontsize}{Node label size in points (default \code{7}).}
#'   }
#'
#' @return An HTML widget created by \code{\link[networkD3]{sankeyNetwork}}. It
#'   renders when printed, in the RStudio viewer or a browser, and can be written
#'   to a standalone file with \code{htmlwidgets::saveWidget}.
#'
#' @import networkD3
#' @importFrom jsonlite toJSON
#'
#' @export
#'
#' @examples
#' # prepare data (a small subset keeps the example quick)
#' G <- sim_data$G[1:150, ]
#' Z <- sim_data$Z[1:150, ]
#' Y_normal <- sim_data$Y_normal[1:150, , drop = FALSE]
#'
#' # plot lucid model
#' fit1 <- estimate_lucid(G = G, Z = Z, Y = Y_normal, lucid_model = "early",
#' CoY = NULL, family = "normal", K = 2, seed = 1008,
#' max_itr = 20, max_tot.itr = 50)
#' plot(fit1)
#'
#' # change node color
#' plot(fit1, G_color = "yellow")
#' plot(fit1, Z_color = "red")
#'
#' # change link color
#' plot(fit1, pos_link_color = "red", neg_link_color = "green")
plot.early_lucid <- function(x, ...) {
  if (!inherits(x, "early_lucid")) {
    stop("`x` must be a fitted early_lucid model (from lucid() or estimate_lucid() ",
        "with lucid_model = 'early').")
  }
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

#' Sankey diagram for a serial-integration LUCID model (not yet implemented)
#'
#' @param x A LUCID model fitted with \code{lucid_model = "serial"}.
#' @param ... Accepted for consistency with \code{\link{plot.early_lucid}}'s
#'   appearance options, but unused.
#' @return Does not return: always raises an error.
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

#' Sankey diagram for a parallel-integration LUCID model (not yet implemented)
#'
#' @param x A LUCID model fitted with \code{lucid_model = "parallel"}.
#' @param ... Accepted for consistency with \code{\link{plot.early_lucid}}'s
#'   appearance options, but unused.
#' @return Does not return: always raises an error.
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

#' @title Network reconstruction from PPCA
#' 
#' @description Constructs a conditional independence network
#'   of the observed variables from the data using the implicitly
#'   estimated covariance matrix within PPCA.
#'
#' @param ppcaOutput \code{list} -- the output object from running any
#'   of the PPCA functions in this package.
#' @param plot \code{logical} -- visualise the resulting network.
#' @param verbose \code{logical} -- verbose intermediary output.
#' @param vertex.size see \code{\link[igraph]{igraph.plotting}}
#' @param edge.width see \code{\link[igraph]{igraph.plotting}}
#' @param vertex.label.cex see \code{\link[igraph]{igraph.plotting}}
#' @param vertex.color see \code{\link[igraph]{igraph.plotting}}
#' @param vertex.label.color see \code{\link[igraph]{igraph.plotting}}
#' @param edge.color see \code{\link[igraph]{igraph.plotting}}
#' @param vertex.label.family see \code{\link[igraph]{igraph.plotting}}
#' @param vertex.label see \code{\link[igraph]{igraph.plotting}}
#' 
#' @details Covariance estimation is done as a preliminary step for
#'   this function. The function then inverts this matrix, which can
#'   be done very efficiently, to obtain the precision matrix. Then
#'   the precision matrix is scaled to unit variance (diagonal) to
#'   obtain partial correlation estimates in the off-diagonal entries,
#'   which is a measure of conditional independence. A two component
#'   mixture model is then fit to the distribution of partial 
#'   correlations using \code{\link[fdrtool:fdrtool]{fdrtool}}. The
#'   partial correlations that are not part of the 'null' component
#'   are then selected as true edges of the network, effectively 
#'   setting the null values to 0. The function then
#'   visualises the resulting network using \code{\link[igraph]{plot.igraph}}.
#'   The user can extract the \code{fdr.stats} element of this output to
#'   view the full output of \code{\link[fdrtool:fdrtool]{fdrtool}},
#'   from which the magnitude and significance of each partial 
#'   correlation can be seen (and customised thresholding can be
#'   performed). The \code{graph} element of the output is an
#'   \linkS4class{igraph} object, and so can be used to easily
#'   make alternative visualisations or compute graph statistics.
#'
#' @return {A \code{list} of 2 elements:
#' \describe{
#' \item{graph}{\linkS4class{igraph} -- Contains the network 
#' information.}
#' \item{fdr.stats}{\code{list} -- the full output of an internal call
#' to \code{\link[fdrtool:fdrtool]{fdrtool}}. Can be useful to inspect
#' the statistics upon which the network was reconstructed.}
#' }}
#' @export
#'
#' @examples
ppca2Net <- function(ppcaOutput, plot=TRUE, verbose=TRUE, vertex.size = 10, 
                     edge.width = 2, vertex.label.cex = 0.4, 
                     vertex.color = "cyan", vertex.label.color = "black", 
                     edge.color = "pink", vertex.label.family= "Helvetica", 
                     vertex.label = NULL){
  
  # extract the covariance matrix estimate
  cov.mat <- ppcaOutput$Sigma
  
  # get variable names
  if(!is.null(colnames(cov.mat))){
    labs <- colnames(cov.mat)
  } else {
    labs <- NULL
    if(verbose){
      cat("no names for variables given... \n")
    }
  }
  
  # extract the partial correlation matrix
  if(verbose){
    cat("computing partial correlations... \n")
  }
  pcor.mat <- -cov2cor(ppca2Covinv(ppcaOutput))
  diag(pcor.mat) <- rep(1, ncol(pcor.mat))
  
  # use fdrtool to assess significance
  if(verbose){
    cat("calculating edge statistics... \n")
  }
  
  edge.stats <- fdrtool::fdrtool(x = c(pcor.mat[lower.tri(pcor.mat)]),
                                 statistic = "correlation",
                                 plot = FALSE,
                                 color.figure = FALSE,
                                 verbose = FALSE)
  
  # create an igraph object
  if(verbose){
    cat("creating network... \n")
  }
  if(!any(edge.stats$qval<edge.stats$param[1])){
    warning("no significant edges at this threshold: inspect $fdr.stats")
  }
  net.thresh <- ifelse(edge.stats$qval<edge.stats$param[1],
                       pcor.mat[lower.tri(pcor.mat)],0)
  net.mat <- matrix(0,
                    nrow = nrow(pcor.mat),
                    ncol = ncol(pcor.mat))
  net.mat[lower.tri(net.mat)] <- net.thresh
  net.mat <- net.mat + t(net.mat)
  diag(net.mat) <- rep(1, nrow(pcor.mat))
  colnames(net.mat) <- labs
  rownames(net.mat) <- labs
  
  
  g <- igraph::graph_from_adjacency_matrix(adjmatrix = net.mat,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           diag = FALSE)
  g <- igraph::delete.vertices(g, igraph::degree(g)==0)
  
  #nEdges     <- gsize(g)
  
  
  igraph::V(g)$color        <- vertex.color #define node colour
  igraph::V(g)$frame.color  <- igraph::V(g)$color #define node colour
  igraph::V(g)$size         <- vertex.size #define node colour
  if(!is.null(vertex.label))
  {
    igraph::V(g)$label        <- vertex.label
  }
  igraph::E(g)$color        <- edge.color#edgeColors
  igraph::E(g)$width        <- edge.width*(abs(igraph::E(g)$weight)/max(abs(igraph::E(g)$weight)))
  igraph::V(g)$label.family <- vertex.label.family
  igraph::V(g)$label.cex    <- vertex.label.cex
  igraph::V(g)$label.color  <- vertex.label.color 
  
  # look to plot the graph
  if (plot && any(edge.stats$qval<edge.stats$param[1])){
    plot(g,
         #main = "Estimated network",
         #vertex.size = vertex.size,
         #edge.width = edge.width*(E(g)$weight/max(E(g)$weight)),
         #vertex.label.family = "Helvetica",
         #vertex.frame.color="red",
         #edge.color = edgeColors,
         #main = "Conditional dependence network",
         layout=igraph::layout.fruchterman.reingold
    )
  }
  
  # create output
  net <- list()
  net[["graph"]] <- g
  net[["fdr.stats"]] <- edge.stats
  
  return(net)
}
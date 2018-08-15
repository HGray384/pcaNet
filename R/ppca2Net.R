ppca2Net <- function(ppcaOutput, plot=TRUE, verbose=TRUE){
  
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
  pcor.mat <- -cov2cor(chol2inv(chol(cov.mat)))
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
  
  # look to plot the graph
  if (plot && any(edge.stats$qval<edge.stats$param[1])){
    plot(g,
         main = "Estimated network",
         vertex.size = 10,
         edge.width = 2)
  }
  
  # create output
  net <- list()
  net[["graph"]] <- g
  net[["fdr.stats"]] <- edge.stats
  
  return(net)
}
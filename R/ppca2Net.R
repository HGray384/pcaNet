ppca2Net <- function(cov.mat, method.ggm=c("prob", "qval","number"),
                     cutoff.ggm=0.8, plot=TRUE, verbose=TRUE){
  
  # get variable names
  if(!is.null(colnames(cov.mat))){
    labs <- colnames(cov.mat)
  } else {
    labs <- 1:ncol(cov.mat)
    if(verbose){
      cat("no names for variables given \n")
    }
  }
  
  # extract the partial correlation matrix
  if(verbose){
    cat("computing partial correlations \n")
  }
  pcor.mat <- corpcor::cor2pcor(cov2cor(cov.mat))
  
  # use GeneNet to assess significance
  if(verbose){
    cat("calculating edge and network statistics \n")
  }
  edges <- GeneNet::network.test.edges(pcor.mat, fdr = TRUE, direct = FALSE,
                                       plot=plot)
  network <- GeneNet::extract.network(edges, 
                                      method.ggm=method.ggm,
                                      cutoff.ggm = cutoff.ggm,
                                      verbose = verbose)
  
  # get graphNEL object
  if(verbose){
    cat("creating graphNEL object \n")
  }
  network.nel <- GeneNet::network.make.graph(network,
                                                  node.labels = labs,
                                                  drop.singles = FALSE)
  
  # # get the adjacency matrices
  # if(verbose){
  #   cat("creating graphAM object (adjacency matrix) \n")
  # }
  # network.AM <- as(network.nel, "graphAM")

  # create output
  net <- list()
  net[["edge.stats"]] <- edges
  net[["net.stats"]] <- network.nel
  
  return(net)
}
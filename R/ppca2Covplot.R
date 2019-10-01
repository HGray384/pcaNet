#' Covariance matrix visualisation
#'
#' @description Heatmap visualisation of the covariance matrix
#'    estimated within PPCA.
#'
#' @param ppcaOutput \code{list} -- the output object from running any
#'   of the PPCA functions in this package.
#'
#' @return plot of the estimated covariance matrix
#'
#' @export
#'
#' @examples
#' #' # simulate a dataset from a zero mean factor model X = Wz + epsilon
#' # start off by generating a random binary connectivity matrix
#' n.factors <- 5
#' n.genes <- 200
#' # with dense connectivity
#' # set.seed(20)
#' conn.mat <- matrix(rbinom(n = n.genes*n.factors,
#'                           size = 1, prob = 0.7), c(n.genes, n.factors))
#' 
#' # now generate a loadings matrix from this connectivity
#' loading.gen <- function(x){
#'   ifelse(x==0, 0, rnorm(1, 0, 1))
#' }
#' 
#' W <- apply(conn.mat, c(1, 2), loading.gen)
#' 
#' # generate factor matrix
#' n.samples <- 100
#' z <- replicate(n.samples, rnorm(n.factors, 0, 1))
#' 
#' # generate a noise matrix
#' sigma.sq <- 0.1
#' epsilon <- replicate(n.samples, rnorm(n.genes, 0, sqrt(sigma.sq)))
#' 
#' # by the ppca equations this gives us the data matrix
#' X <- W%*%z + epsilon
#' WWt <- tcrossprod(W)
#' Sigma <- WWt + diag(sigma.sq, n.genes)
#' 
#' # select 10% of entries to make missing values
#' missFrac <- 0.1
#' inds <- sample(x = 1:length(X),
#'                size = ceiling(length(X)*missFrac),
#'                replace = FALSE)
#' 
#' # replace them with NAs in the dataset
#' missing.dataset <- X
#' missing.dataset[inds] <- NA
#' 
#' # run ppca
#' ppf <- pca_full(missing.dataset, ncomp=5, algorithm="vb", maxiters=5,
#' bias=TRUE, rotate2pca=FALSE, loglike=TRUE, verbose=TRUE)
#' 
#' # plot the matrix
#' ppca2Covplot(ppf)
ppca2Covplot <- function(ppcaOutput){
  paletteLength <- 50
  myColor <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                              "PuOr")))(paletteLength)
  myMax <- max(max(ppcaOutput$Sigma), abs(min(ppcaOutput$Sigma)))
  myBreaks <- c(seq(-myMax, 0, length.out=ceiling(paletteLength/2) + 1),
                seq(myMax/paletteLength, max(ppcaOutput$Sigma),
                    length.out=floor(paletteLength/2)))
  pheatmap::pheatmap(ppcaOutput$Sigma, color=myColor, breaks=myBreaks, fontsize =
                       8, cluster_rows = F, cluster_cols = F, legend = T, show_rownames = F,
                     show_colnames = F, border_color = NA)
}
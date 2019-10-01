#' Inverse covariance matrix computation from PPCA
#'
#' @description Efficient inversion of the covariance matrix
#'   estimated from PPCA. 
#'   
#' @param ppcaOutput \code{list} -- the output object from running any
#'   of the PPCA functions in this package.
#'
#' @details The computation exploits the
#'   Woodbury identity so that a kxk matrix (where k is often less than 10)
#'   is inverted instead of the potentially large pxp matrix. The closed-form
#'   expression for the inverse depends upon parameters that are estimated
#'   in the PPCA algorithm.
#'
#' @return \code{matrix} -- the inverse of the covariance matrix.
#' @export
#'
#' @examples
#' # simulate a dataset from a zero mean factor model X = Wz + epsilon
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
#' # compute the inverse
#' covinv <- ppca2Covinv(ppf)
#' system.time(ppca2Covinv(ppf))
#' 
#' covinv2 <- solve(ppf$Sigma)
#' system.time(solve(ppf$Sigma))
ppca2Covinv <- function(ppcaOutput){
  # get the relevant parameters from the ppca output
  sigmaSq <- ppcaOutput$sigmaSq
  W <- ppcaOutput$W
  k <- ncol(W)
  p <- nrow(W)
  
  # compute the kxk matrix to be inverted and invert it
  M <- sigmaSq*diag(k)+crossprod(W)
  Minv <- chol2inv(chol(M))
  
  # compute the covariance inverse and return
  # the formula for the inversion comes from applying the
  # Woodbury formula to the covariance matrix decomposition
  # of ppca
  sigmaInv <-  (1/sigmaSq)*(diag(p)-W%*%tcrossprod(Minv,W))
  return(sigmaInv)
}
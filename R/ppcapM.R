#' Implements a probabilistic PCA missing value estimator, as in pcaMethods.
#'   Use of Rcpp makes this version faster (cite software note) and 
#'   the emphasised output is the covariance matrix \code{Sigma}, which 
#'   can be used for network reconstruction.
#' 
#' Details about the probabilistic model underlying PPCA are found in
#' Bishop 1999. The algorithm (Verbeek, 20??) uses an expectation maximation
#' approach together with a probabilistic model to approximate the
#' principal axes (eigenvectors of the covariance matrix in PCA).
#' The estimation is done iteratively, the algorithm terminates if
#' either the maximum number of iterations is reached or if the
#' estimated increase in precision falls below \eqn{1e^{-4}}{1e^-4}.
#' 
#' @title Probabilistic PCA (pcaMethods version)
#'
#' @param myMat \code{matrix} -- Pre-processed matrix (centered,
#'   scaled) with variables in columns and observations in rows. The
#'   data may contain missing values, denoted as \code{NA}. 
#' @param nPcs \code{numeric} -- Number of components used for
#'   re-estimation. Choosing few components may decrease the
#'   estimation precision.
#' @param seed \code{numeric} -- the random number seed used, useful
#'   to specify when comparing algorithms.
#' @param threshold \code{numeric} -- Convergence threshold. 
#'   If the increase in precision of an update
#'   falls below this then the algorithm is stopped.
#' @param maxIterations  \code{numeric} -- Maximum number of estimation
#'   steps. 
#' @param ... 
#'
#' @return {A \code{list} of 4 elements:
#' \describe{
#' \item{W}{\code{matrix} -- the estimated loadings.}
#' \item{sigmaSq}{\code{numeric} -- the estimated isotropic variance.}
#' \item{Sigma}{\code{matrix} -- the estimated covariance matrix.}
#' \item{pcaMethodsRes}{\code{class} -- 
#'   see \code{pcaRes}.}
#' }}
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
#'                replace = F)
#' 
#' # replace them with NAs in the dataset
#' missing.dataset <- X
#' missing.dataset[inds] <- NA
#' 
#' # run ppca
#' pp <- ppcapM(t(missing.dataset), nPcs = 5)
#' names(pp)
#' 
#' # sigmasq estimation
#' abs(pp$sigmaSq-sigma.sq)
#' 
#' # X reconstruction
#' recon.X <- pp$pcaMethodsRes@loadings%*%t(pp$pcaMethodsRes@scores)
#' norm(recon.X-X, type="F")^2/(length(X))
#' 
#' # covariance estimation
#' norm(pp$Sigma-Sigma, type="F")^2/(length(X))
ppcapM <- function(myMat, nPcs=2, seed=NA, threshold=1e-4, maxIterations=1000, ...) {

  if (!is.na(seed)) 
    set.seed(seed)
  
  N <- nrow(myMat)
  D <- ncol(myMat)
  
  hidden   <- which(is.na(myMat))
  nMissing <- length(hidden)
  
  myMatsaved      <- myMat
  
  
  if(nMissing) { myMat[hidden] <- 0 } 
  
  ## ------- Initialization
  r <- sample(N)
  C <- t(myMat[r[1:nPcs], ,drop = FALSE])
  ## Random matrix with the same dimnames as myMat
  C <- matrix(rnorm(C), nrow(C), ncol(C), dimnames = labels(C) )
  
  ppcaOutput <- ppcaNet(myMat, N, D, C,  hidden, nMissing, nPcs, threshold, maxIterations=1000)
  
  myMat <- ppcaOutput$myMat;
  
  # Additional processing from pcaMethods to orthonormalise:
  Worth <- pcaMethods:::orth(ppcaOutput$W)
  evs   <- eigen(cov(myMat %*% Worth))
  vals  <- evs[[1]]
  vecs  <- evs[[2]]
  Worth <- Worth %*% vecs
  X     <- myMat %*% Worth
  R2cum <- rep(NA, nPcs)
  TSS   <- sum(myMat^2, na.rm = TRUE)
  for (i in 1:ncol(Worth)) {
    difference <- myMat - (X[, 1:i, drop = FALSE] %*% t(Worth[, 1:i, drop = FALSE]))
    R2cum[i] <- 1 - (sum(difference^2, na.rm = TRUE)/TSS)
  }
  
  
  # Prepare pcaMethods-style output:
  pcaMethodsRes          <- new("pcaRes")
  pcaMethodsRes@scores   <- X
  pcaMethodsRes@loadings <- Worth
  pcaMethodsRes@R2cum    <- R2cum
  pcaMethodsRes@method   <- "ppca"
  
  # Also return the standard ppcaNet output:
  output <- list()
  output[["W"]]              <- ppcaOutput$W
  output[["sigmaSq"]]        <- ppcaOutput$ss
  output[["Sigma"]]          <- ppcaOutput$C
  output[["pcaMethodsRes"]]  <- pcaMethodsRes
  
  return(output)
  
}
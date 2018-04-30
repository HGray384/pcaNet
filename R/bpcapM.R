#' Implements a Bayesian PCA missing value estimator, as in pcaMethods.
#' 
#' Details about the probabilistic model underlying BPCA are found in
#' Oba et. al 2003. The algorithm uses an expectation maximation
#' approach together with a Bayesian model to approximate the
#' principal axes (eigenvectors of the covariance matrix in PCA).
#' The estimation is done iteratively, the algorithm terminates if
#' either the maximum number of iterations was reached or if the
#' estimated increase in precision falls below \eqn{1e^{-4}}{1e^-4}.
#' 
#' @title Bayesian PCA missing value estimation (pcaMethods version)
#'
#' @param myMat \code{matrix} -- Pre-processed matrix (centered,
#'   scaled) with variables in columns and observations in rows. The
#'   data may contain missing values, denoted as \code{NA}. 
#' @param nPcs \code{numeric} -- Number of components used for
#'   re-estimation. Choosing few components may decrease the
#'   estimation precision.
#' @param threshold convergence threshold.
#' @param maxIterations  \code{numeric} -- Maximum number of estimation
#'   steps. 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
bpcapM <- function(myMat, nPcs=2, threshold=1e-4, maxIterations=100, ...) {

  N <- nrow(myMat)
  D <- ncol(myMat)
  
  if (!is.na(nPcs)) 
  {
    qDim <- nPcs
  }
  else
  {
    qDim <- D - 1
  }
  
  
  isNAmat <- is.na(myMat)
  numberOfNonNAvaluesInEachCol <- colSums(!isNAmat)
  numberOfNonNAvaluesInEachRow <- rowSums(!isNAmat)

  hidden          <- which(isNAmat)
  nomissIndex     <- which(numberOfNonNAvaluesInEachRow == D)
  missIndex       <- which(numberOfNonNAvaluesInEachRow != D)
  
  myMatsaved      <- myMat
  
  nMissing <- length(hidden)
  if(nMissing) { 
    myMat[hidden] <- 0 
  } 
  
  print(dim(myMat))
  print(hidden)
  
  ppcaOutput <- bpcaNet(myMat, N, D, hidden, numberOfNonNAvaluesInEachCol, nomissIndex, missIndex, nMissing, qDim, threshold, maxIterations=1000)

  
  R2cum      <- rep(NA, nPcs)
  TSS        <- sum(myMatsaved^2, na.rm=TRUE)
  
  for (i in 1:nPcs) {
    difference <- myMatsaved - (ppcaOutput$scores[,1:i, drop=FALSE] %*% t(ppcaOutput$W[,1:i, drop=FALSE]) )
    R2cum[i]   <- 1 - (sum(difference^2, na.rm=TRUE) / TSS)
  }
  
  
  pcaMethodsRes           <- new("pcaRes")
  pcaMethodsRes@scores    <- ppcaOutput$scores 
  pcaMethodsRes@loadings  <- ppcaOutput$W
  pcaMethodsRes@R2cum     <- R2cum
  pcaMethodsRes@method    <- "bpca"
  
  
  # Return standard ppcaNet output:
  
  output <- list()
  output[["W"]]              <- ppcaOutput$W
  output[["sigmaSq"]]        <- ppcaOutput$ss
  output[["Sigma"]]          <- ppcaOutput$C
  output[["pcaMethodsRes"]]  <- pcaMethodsRes
  

  return(output)
}
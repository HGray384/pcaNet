#' Title
#'
#' @param myMat 
#' @param nPcs 
#' @param maxIterations 
#' @param TolFun 
#' @param TolX 
#'
#' @return
#' @export
#'
#' @examples
mpca <- function(myMat, nPcs=2, maxIterations = 10000 ,TolFun = 1e-6,TolX = 1e-6) {

  k        <- nPcs
  myMat    <- t(myMat)
  p        <- nrow(myMat)
  n        <- ncol(myMat)
  
  traceS   <- sum(myMat^2)/(n-1)
  
  hidden   <- which(is.na(myMat))
  nMissing <- length(hidden)
  
  if(nMissing) { myMat[hidden] <- 0 } 
  
  ## Initialization
  W        <- matrix(rnorm(p*k), nrow = p, ncol = k)
  v        <- runif(1)

  ## Call cpp code to perform processing:
  ppcaOutput  <- mpcaNet(myMat, W, hidden, nMissing, v, traceS, maxIterations, TolFun, TolX)
  
  if(ppcaOutput$numIter == maxIterations)
  {
    print('Maximum number of iterations reached')
  }
  
  output <- list()
  output[["W"]]        <- ppcaOutput$W
  output[["sigmaSq"]]  <- ppcaOutput$ss
  output[["Sigma"]]    <- ppcaOutput$C
  output[["numIter"]]  <- ppcaOutput$numIter
  
  return(output)
  
}
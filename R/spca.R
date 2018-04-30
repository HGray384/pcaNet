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
spcapp <- function (myMat, nPcs, maxIterations = 10000 ,TolFun = 1e-6,TolX = 1e-6)
{
    k        <- nPcs
    myMat    <- scale(myMat, scale = F)
    myMat    <- t(myMat)
    p        <- nrow(myMat)
    n        <- ncol(myMat)
    
    traceS   <- sum(myMat^2)/(n-1)

    # Initialise W and v    
    W        <- matrix(rnorm(p*k), nrow = p, ncol = k)
    v        <- runif(1)
    
    ppcaOutput  <- ppcaSensible(myMat, W, v, traceS, maxIterations, TolFun, TolX)
    
    if(ppcaOutput$numIter == maxIterations)
    {
        print('Maximum number of iterations reached')
    }
    
    # Return standard ppcaNet output:
    
    output <- list()
    output[["W"]]        <- ppcaOutput$W
    output[["sigmaSq"]]  <- ppcaOutput$ss
    output[["Sigma"]]    <- ppcaOutput$C
    output[["numIter"]]  <- ppcaOutput$numIter
    
    return(output)
}

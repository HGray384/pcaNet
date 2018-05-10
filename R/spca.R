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
spca <- function (myMat, nPcs = 2, maxIterations = 10000 ,TolFun = 1e-4,TolX = 1e-4)
{
    k        <- nPcs
    myMat    <- scale(myMat, scale = F)
    myMatsaved <- myMat
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
    
    R2cum      <- rep(NA, nPcs)
    TSS        <- sum(myMatsaved^2)
    
    for (i in 1:nPcs) {
      difference <- myMatsaved - (ppcaOutput$scores[,1:i, drop=FALSE] %*% t(ppcaOutput$W[,1:i, drop=FALSE]) )
      R2cum[i]   <- 1 - (sum(difference^2) / TSS)
    }
    
    pcaMethodsRes           <- new("pcaRes")
    pcaMethodsRes@scores    <- ppcaOutput$scores 
    pcaMethodsRes@loadings  <- ppcaOutput$W
    pcaMethodsRes@R2cum     <- R2cum
    pcaMethodsRes@method    <- "spca"
    
    # Return standard ppcaNet output:
    
    output <- list()
    output[["W"]]              <- ppcaOutput$W
    output[["sigmaSq"]]        <- ppcaOutput$ss
    output[["Sigma"]]          <- ppcaOutput$C
    output[["numIter"]]  <- ppcaOutput$numIter
    output[["pcaMethodsRes"]]  <- pcaMethodsRes
    
    return(output)
}

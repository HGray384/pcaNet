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
vbpca <- function(myMat, nPcs=NA, maxIterations = 10000 ,TolFun = 1e-4,TolX = 1e-4) {
  
  p        <- nrow(myMat)
  n        <- ncol(myMat)
  if(is.na(nPcs)){
    nPcs <- n-1
  }
  k        <- nPcs
  myMatsaved <- myMat


  
  traceS   <- sum(myMat^2)/(n-1)
  
  hidden   <- which(is.na(myMat))
  nMissing <- length(hidden)
  
  if(nMissing) { myMat[hidden] <- 0 } 
  
  ## ------- Initialization
  #r <- sample(n)
  #C <- t(myMat[r[1:nPcs], ,drop = FALSE])
  ## Random matrix with the same dimnames as myMat
  #C <- matrix(rnorm(C), nrow(C), ncol(C), dimnames = labels(C) )
  #W <- matrix(c(-1.900198139320733, 0.683830133135421, -0.552842656796558,-0.672586670748523,-0.536357706683455,-0.577405592794564,-0.994303992040422,-0.463725494945621,-0.362414139913114,0.158209008970325,2.441328431935101,-1.385253442492738), nrow = 4, ncol = 3, byrow = T)
  #v <- 0.048375054285105
  #ppcaOutput <- mpcaNet(myMat, N, D, C,  hidden, nMissing, nPcs, threshold, maxIterations=1000)
  #myMat <- t(myMat)
  W        <- qr.Q(qr(matrix(rnorm(p*k), nrow = p, ncol = k))) # need orthonormal basis so covariances=0
  v        <- runif(1)
  
  ppcaOutput  <- vbpcaNet(myMat, W, hidden, nMissing, v, traceS, maxIterations, TolFun, TolX)
  
  if(ppcaOutput$numIter == maxIterations)
  {
    print('Maximum number of iterations reached')
  }
  R2cum      <- rep(NA, nPcs)
  TSS        <- sum(myMatsaved^2, na.rm = TRUE)
  
  for (i in 1:nPcs) {
    difference <- t(myMatsaved) - (ppcaOutput$scores[,1:i, drop=FALSE] %*% t(ppcaOutput$W[,1:i, drop=FALSE]) )
    R2cum[i]   <- 1 - (sum(difference^2, na.rm = TRUE) / TSS)
  }
  
  pcaMethodsRes           <- new("pcaRes")
  pcaMethodsRes@scores    <- ppcaOutput$scores 
  pcaMethodsRes@loadings  <- ppcaOutput$W
  pcaMethodsRes@R2cum     <- R2cum
  pcaMethodsRes@method    <- "vbpca"
  
  # Return standard ppcaNet output:
  
  output <- list()
  output[["W"]]              <- ppcaOutput$W
  output[["sigmaSq"]]        <- ppcaOutput$ss
  output[["Sigma"]]          <- ppcaOutput$C
  output[["numIter"]]  <- ppcaOutput$numIter
  output[["pcaMethodsRes"]]  <- pcaMethodsRes
  
  return(output)
}
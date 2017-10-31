vbpcapp <- function(Y, nPcs=2, MaxIter = 10000 ,TolFun = 1e-6,TolX = 1e-6) {
  
  k        <- nPcs
  #Y        <- scale(Y, scale = F)
  Y        <- t(Y)
  p        <- nrow(Y)
  n        <- ncol(Y)
  
  traceS   <- sum(Y^2)/(n-1)
  
  hidden   <- which(is.na(Y))
  nMissing <- length(hidden)
  
  if(nMissing) { Y[hidden] <- 0 } 
  
  ## ------- Initialization
  #r <- sample(n)
  #C <- t(myMat[r[1:nPcs], ,drop = FALSE])
  ## Random matrix with the same dimnames as myMat
  #C <- matrix(rnorm(C), nrow(C), ncol(C), dimnames = labels(C) )
  #W <- matrix(c(-1.900198139320733, 0.683830133135421, -0.552842656796558,-0.672586670748523,-0.536357706683455,-0.577405592794564,-0.994303992040422,-0.463725494945621,-0.362414139913114,0.158209008970325,2.441328431935101,-1.385253442492738), nrow = 4, ncol = 3, byrow = T)
  #v <- 0.048375054285105
  #ppcaOutput <- mpcaNet(myMat, N, D, C,  hidden, nMissing, nPcs, threshold, maxIterations=1000)
  #Y <- t(myMat)
  W        <- qr.Q(qr(matrix(rnorm(p*k), nrow = p, ncol = k))) # need orthonormal basis so covariances=0
  v        <- runif(1)
  
  ppcaOutput  <- vbpcaNet(Y, W, hidden, nMissing, v, traceS, MaxIter, TolFun, TolX)
  
  if(ppcaOutput$numIter == MaxIter)
  {
    print('Maximum number of iterations reached')
  }
  res <- list()
  res[["W"]]        <- ppcaOutput$W
  res[["sigmaSq"]]  <- ppcaOutput$ss
  res[["C"]]        <- ppcaOutput$C
  res[["method"]]   <- "ppcaNet"
  return(res)
  
}
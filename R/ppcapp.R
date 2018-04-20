ppcapp <- function(myMat, nPcs=2, seed=NA, threshold=1e-5, maxIterations=1000, ...) {
  ## Set the seed to the user defined value. This affects the generation
  ## of random values for the initial setup of the loading matrix
  
  if (!is.na(seed)) 
    set.seed(seed)
  
  N <- nrow(myMat)
  D <- ncol(myMat)
  
  hidden   <- which(is.na(myMat))
  nMissing <- length(hidden)
  
  if(nMissing) { myMat[hidden] <- 0 } 
  
  ## ------- Initialization
  r <- sample(N)
  C <- t(myMat[r[1:nPcs], ,drop = FALSE])
  ## Random matrix with the same dimnames as myMat
  C <- matrix(rnorm(C), nrow(C), ncol(C), dimnames = labels(C) )
  
  ppcaOutput <- ppcaNet(myMat, N, D, C,  hidden, nMissing, nPcs, threshold, maxIterations=1000)
  
  Worth <- pcaMethods::orth(ppcaOutput$W)
  evs   <- eigen(cov(myMat %*% C))
  vals  <- evs[[1]]
  vecs  <- evs[[2]]
  Worth <- Worth %*% vecs
  X     <- myMat %*% Worth
  R2cum <- rep(NA, nPcs)
  TSS   <- sum(myMat^2, na.rm = TRUE)
  for (i in 1:ncol(Worth)) {
    difference <- myMat - (X[, 1:i, drop = FALSE] %*% t(C[, 1:i, drop = FALSE]))
    R2cum[i] <- 1 - (sum(difference^2, na.rm = TRUE)/TSS)
  }
  
  
  
  res <- list()
  res[["W"]]        <- ppcaOutput$W
  res[["W_orth"]]   <- Worth
  res[["sigmaSq"]]  <- ppcaOutput$ss
  res[["C"]]        <- ppcaOutput$C
  res[["method"]]   <- "ppcaNet"
  return(res)
  
}
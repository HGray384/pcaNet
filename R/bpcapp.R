bpcapp <- function(myMat, nPcs=NA, seed=NA, threshold=1e-5, maxIterations=1000, ...) {
  ## Set the seed to the user defined value. This affects the generation
  ## of random values for the initial setup of the loading matrix
  
  if (!is.na(seed)) 
    set.seed(seed)
  

  
  
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
  
  print(qDim)
  
  isNAmat <- is.na(myMat)
  numberOfNonNAvaluesInEachCol <- colSums(!isNAmat)
  numberOfNonNAvaluesInEachRow <- rowSums(!isNAmat)
  print(numberOfNonNAvaluesInEachCol)
  print(numberOfNonNAvaluesInEachRow)
  hidden          <- which(isNAmat)
  nomissIndex     <- which(numberOfNonNAvaluesInEachRow == D)
  missIndex       <- which(numberOfNonNAvaluesInEachRow != D)
  print(missIndex)
  #hiddenArrayInds <- arrayInd(which(is.na(X)), dim(X))
  nMissing <- length(hidden)
  if(nMissing) { 
    myMat[hidden] <- 0 
  } 
  
  t <- bpcaNet(myMat, N, D, hidden, numberOfNonNAvaluesInEachCol, nomissIndex, missIndex, nMissing, qDim, threshold, maxIterations=1000)

  res <- list()
  res[["W"]]        <- t$C
  res[["sigmaSq"]]  <- t$ss
  res[["scores"]]   <- t$scores
  res[["loadings"]] <- t$loadings
  res[["R2cum"]]    <- t$R2cum
  res[["method"]]   <- "ppca"
  return(res)
  
}
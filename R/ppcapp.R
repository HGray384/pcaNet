ppcapp <- function(myMat, nPcs=2, seed=NA, threshold=1e-5, maxIterations=1000, ...) {
  ## Set the seed to the user defined value. This affects the generation
  ## of random values for the initial setup of the loading matrix
  
  if (!is.na(seed)) 
    set.seed(seed)
  
  N <- nrow(myMat)
  D <- ncol(myMat)
  
  hidden <- which(is.na(myMat))
  nMissing <- length(hidden)
  
  if(nMissing) { myMat[hidden] <- 0 } 
  
  ## ------- Initialization
  r <- sample(N)
  C <- t(myMat[r[1:nPcs], ,drop = FALSE])
  ## Random matrix with the same dimnames as myMat
  C <- matrix(rnorm(C), nrow(C), ncol(C), dimnames = labels(C) )
  
  t <- ppcaNet(myMat, N, D, C,  hidden, nMissing, nPcs, threshold, maxIterations=1000)

  res <- list()
  res[["W"]]        <- t$C
  res[["sigmaSq"]]  <- t$ss
  res[["scores"]]   <- t$scores
  res[["loadings"]] <- t$loadings
  res[["R2cum"]]    <- t$R2cum
  res[["method"]]   <- "ppca"
  return(res)
  
}
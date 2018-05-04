#' Title
#'
#' @param myMat 
#' @param nPcs 
#' @param seed 
#' @param threshold 
#' @param maxIterations 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
ppcapM <- function(myMat, nPcs=2, seed=NA, threshold=1e-5, maxIterations=1000, ...) {

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
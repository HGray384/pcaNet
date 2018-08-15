pcapM <- function(myMat, nPcs=2, method='ppca', seed=NA, threshold=1e-4,
                  maxIterations=1000, center = TRUE, scale = TRUE, 
                  loglike = TRUE, verbose=TRUE, ...) {
  # preprocessing
  if (nPcs > ncol(myMat)) {
    warning("more components than matrix columns requested")
    nPcs <- min(dim(myMat))
  }
  if (nPcs > nrow(myMat)) {
    warning("more components than matrix rows requested")
    nPcs <- min(dim(myMat))
  }
  
  missing <- is.na(myMat)
  myMat <- scale(myMat, center = center, scale = scale)
  m <- attr(myMat, "scaled:center")
  sc <- attr(myMat, "scaled:scale")
  if (is.null(m)){
    m <- rep(0, ncol(myMat))
  }
  if (is.null(sc)){
    sc <- rep(1, ncol(myMat))
  }
  if (scale){
    if (center) {
      scale.meth <- "uv"
    }
    else {
      scale.meth <- "rms"
    }
  }
  
  # call to ppcapM or bpcapM
  if (method=="ppca"){
    res <- ppcapM(myMat, nPcs=nPcs, seed=seed, threshold=threshold,
                  maxIterations=maxIterations, loglike = loglike, 
                  verbose=verbose, ...)
  } else if (method=="bpca"){
    res <- bpcapM(myMat, nPcs=nPcs, threshold=threshold,
                  maxIterations=maxIterations, loglike = loglike, 
                  verbose=verbose, ...)
  } else {
    stop("The specified method must be either 'ppca' or 'bpca'")
  }
  
  # structure output
  res$pcaMethodsRes@nPcs <- nPcs # do we need to edit this?
  res$pcaMethodsRes@nObs <- nrow(myMat)
  res$pcaMethodsRes@nVar <- ncol(myMat)
  res$pcaMethodsRes@sDev <- apply(res$pcaMethodsRes@scores, 2, sd)
  res$pcaMethodsRes@center <- m
  res$pcaMethodsRes@centered <- center
  res$pcaMethodsRes@scale <- sc
  res$pcaMethodsRes@scaled <- scale.meth
  res$pcaMethodsRes@R2 <- res$pcaMethodsRes@R2cum[1]
  if (length(res$pcaMethodsRes@R2cum) > 1) {
    res$pcaMethodsRes@R2 <- c(res$pcaMethodsRes@R2,
                              diff(res$pcaMethodsRes@R2cum))
  }
  
  if(any(missing)){
    completeObs <- myMat
    recData <- tcrossprod(res$pcaMethodsRes@scores[, 1:nPcs, drop = FALSE],
                          res$pcaMethodsRes@loadings[, 1:nPcs, drop = FALSE])
    recData <- sweep(recData, 2, sc, "*")
    recData <- sweep(recData, 2, m, "+")
    completeObs[missing] <- recData[missing]
    res$pcamethodsRes@completeObs <- completeObs
  }
  
  # return results
  return(res)
}
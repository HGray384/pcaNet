#' @title A wrapper for pcaMethods function implementations
#' 
#' @description Implements the equivalent of 
#'   \code{\link[pcaMethods:pca]{pca}}.
#'   This function preprocesses the data as specified by the user,
#'   then calls ppcapM or bpcapM, and finally handles this output
#'   to return a list. One element of the output is a pcaRes object.
#'
#' @param myMat \code{matrix} -- Data matrix with 
#'   variables in columns and observations in rows. The
#'   data may contain missing values, denoted as \code{NA}. 
#' @param nPcs \code{numeric} -- Number of components used for
#'   re-estimation. Choosing few components may decrease the
#'   estimation precision.
#' @param method \code{c("ppca", "bpca")} -- frequentist or
#'   Bayesian estimation of model parameters.
#' @param seed \code{numeric} -- the random number seed used, useful
#'   to specify when comparing algorithms.
#' @param threshold \code{numeric} -- Convergence threshold. 
#'   If the increase in precision of an update
#'   falls below this then the algorithm is stopped.
#' @param maxIterations \code{numeric} -- Maximum number of estimation
#'   steps. 
#' @param center \code{logical} -- should the data be centered?
#' @param scale \code{c("none", "pareto", "vector", "uv")} --
#'   which method of scaling should be used? See 
#'   \code{\link[pcaMethods:pca]{pca}}. Defaults to "none".
#' @param loglike \code{logical} -- should the log-likelihood
#'   of the estimated parameters be returned? See Details.
#' @param verbose \code{logical} -- verbose intermediary 
#'   algorithm output.
#' 
#' @details See \code{\link{ppcapM}} and \code{\link{bpcapM}} for 
#'   the algorithm specifics. \code{loglike} indicates whether 
#'   log-likelihood values for the resulting estimates should 
#'   be computed. This can be useful to compare different algorithms.
#'
#' @return {A \code{list} of 5 or 7 elements, depending on the value
#' of \code{loglike}:
#' \describe{
#' \item{W}{\code{matrix} -- the estimated loadings.}
#' \item{sigmaSq}{\code{numeric} -- the estimated isotropic variance.}
#' \item{Sigma}{\code{matrix} -- the estimated covariance matrix.}
#' \item{m}{\code{numeric} -- the estimated mean vector.}
#' \item{logLikeObs}{\code{numeric} -- the log-likelihood value
#' of the observed data given the estimated parameters.}
#' \item{logLikeImp}{\code{numeric} -- the log-likelihood value
#' of the imputed data given the estimated parameters.}
#' \item{pcaMethodsRes}{\code{class} -- 
#'   see \linkS4class{pcaRes}.}
#' }}
#' @export
#'
#' @examples
#' # simulate a dataset from a zero mean factor model X = Wz + epsilon
#' # start off by generating a random binary connectivity matrix
#' n.factors <- 5
#' n.genes <- 200
#' # with dense connectivity
#' # set.seed(20)
#' conn.mat <- matrix(rbinom(n = n.genes*n.factors,
#'                           size = 1, prob = 0.7), c(n.genes, n.factors))
#' 
#' # now generate a loadings matrix from this connectivity
#' loading.gen <- function(x){
#'   ifelse(x==0, 0, rnorm(1, 0, 1))
#' }
#' 
#' W <- apply(conn.mat, c(1, 2), loading.gen)
#' 
#' # generate factor matrix
#' n.samples <- 100
#' z <- replicate(n.samples, rnorm(n.factors, 0, 1))
#' 
#' # generate a noise matrix
#' sigma.sq <- 0.1
#' epsilon <- replicate(n.samples, rnorm(n.genes, 0, sqrt(sigma.sq)))
#' 
#' # by the ppca equations this gives us the data matrix
#' X <- W%*%z + epsilon
#' WWt <- tcrossprod(W)
#' Sigma <- WWt + diag(sigma.sq, n.genes)
#' 
#' # select 10% of entries to make missing values
#' missFrac <- 0.1
#' inds <- sample(x = 1:length(X),
#'                size = ceiling(length(X)*missFrac),
#'                replace = FALSE)
#' 
#' # replace them with NAs in the dataset
#' missing.dataset <- X
#' missing.dataset[inds] <- NA
#' 
#' # run ppca
#' ppm <- pcapM(t(missing.dataset), nPcs=5, method="bpca", seed=2009, 
#' maxIterations=1000, center=TRUE, loglike=TRUE, verbose=TRUE)
#' 
pcapM <- function(myMat, nPcs=2, method='ppca', seed=NA, threshold=1e-4,
                  maxIterations=1000, center = TRUE, 
                  scale = c("none"), 
                  loglike = TRUE, verbose=TRUE) {
  ## preprocessing
  if (nPcs > ncol(myMat)) {
    warning("more components than matrix columns requested")
    nPcs <- min(dim(myMat))
  }
  if (nPcs > nrow(myMat)) {
    warning("more components than matrix rows requested")
    nPcs <- min(dim(myMat))
  }
  # any missing
  missing <- is.na(myMat)
  
  # scaling and centering
  if(is.null(scale)){
    scale = "none"
  }
  if (length(scale)!=1){
    stop("scale must have length 1")
  } else {
    if(!(scale %in% c("none", "pareto", "vector", "uv"))){
      stop("please provide a valid scaling method")
    }
  }

  if(is.null(center)){
    center = FALSE
  }
  if(!is.logical(center)){
    stop("please provide TRUE or FALSE for centering")
  }
  if(center){
    m <- colMeans(myMat, na.rm = TRUE)
  }
  else {
    m <- rep(0, ncol(myMat))
  }
  myMat <- sweep(myMat, 2, m, "-")
  
  if(scale=="none"){
    sc <- rep(1, ncol(myMat))
  }
  if(scale=="uv"){
    sc <- apply(myMat, 2, sd, na.rm = TRUE)
  }
  if (scale == "pareto") {
    sc <- sqrt(apply(myMat, 2, sd, na.rm = TRUE))
  }
  if (scale == "vector") {
    sc <- apply(myMat, 2, function(x){sqrt(sum(x^2, na.rm = TRUE))})
  }
  myMat <- sweep(myMat, 2, sc, "/")
  
  
  # call to ppcapM or bpcapM
  if (method=="ppca"){
    res <- ppcapM(myMat, nPcs=nPcs, seed=seed, threshold=threshold,
                  maxIterations=maxIterations, loglike = loglike, 
                  verbose=verbose)
  } else if (method=="bpca"){
    res <- bpcapM(myMat, nPcs=nPcs, threshold=threshold,
                  maxIterations=maxIterations, loglike = loglike, 
                  verbose=verbose)
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
  res$pcaMethodsRes@scaled <- scale
  res$pcaMethodsRes@R2 <- res$pcaMethodsRes@R2cum[1]
  if (length(res$pcaMethodsRes@R2cum) > 1) {
    res$pcaMethodsRes@R2 <- c(res$pcaMethodsRes@R2,
                              diff(res$pcaMethodsRes@R2cum))
  }
  
  completeObs <- myMat
  if(any(missing)){
    recData <- tcrossprod(res$pcaMethodsRes@scores[, 1:nPcs, drop = FALSE],
                          res$pcaMethodsRes@loadings[, 1:nPcs, drop = FALSE])
    recData <- sweep(recData, 2, sc, "*")
    recData <- sweep(recData, 2, m, "+")
    completeObs[missing] <- recData[missing]
  }
  res$pcaMethodsRes@completeObs <- completeObs
  
  # return results
  return(res)
}
#' Implements a probabilistic PCA missing value estimator, as in pcaMethods.
#'   Use of Rcpp makes this version faster and 
#'   the emphasised output is the covariance matrix \code{Sigma}, which 
#'   can be used for network reconstruction.
#' 
#' Details about the probabilistic model underlying PPCA are found in
#' Bishop 1999. The algorithm (Porta, 2005) uses an expectation maximisation
#' approach together with a probabilistic model to approximate the
#' principal axes (eigenvectors of the covariance matrix in PCA).
#' The estimation is done iteratively, the algorithm terminates if
#' either the maximum number of iterations is reached or if the
#' estimated increase in precision falls below \eqn{1e^{-4}}{1e^-4}.
#' 
#' @title Probabilistic PCA (pcaMethods version)
#'
#' @param myMat \code{matrix} -- Pre-processed matrix (centered,
#'   scaled) with variables in columns and observations in rows. The
#'   data may contain missing values, denoted as \code{NA}. 
#' @param nPcs \code{numeric} -- Number of components used for
#'   re-estimation. Choosing few components may decrease the
#'   estimation precision.
#' @param seed \code{numeric} -- the random number seed used, useful
#'   to specify when comparing algorithms.
#' @param threshold \code{numeric} -- Convergence threshold. 
#'   If the increase in precision of an update
#'   falls below this then the algorithm is stopped.
#' @param maxIterations  \code{numeric} -- Maximum number of estimation
#'   steps. 
#' @param loglike \code{logical} -- should the log-likelihood
#'   of the estimated parameters be returned? See Details.
#' @param verbose \code{logical} -- verbose intermediary 
#'   algorithm output.
#'
#' @return {A \code{list} of 4 elements:
#' \describe{
#' \item{W}{\code{matrix} -- the estimated loadings.}
#' \item{sigmaSq}{\code{numeric} -- the estimated isotropic variance.}
#' \item{Sigma}{\code{matrix} -- the estimated covariance matrix.}
#' \item{pcaMethodsRes}{\code{class} -- 
#'   see \code{\link[pcaMethods:pcaRes-class]{pcaRes}}.}
#' }}
#' @export
#'
#' @references Porta, J.M., Verbeek, J.J. and Kroese, B.J., 2005.
#'  \href{https://hal.inria.fr/inria-00321476/en}{link}
#' 
#'  Stacklies, W., Redestig, H., Scholz, M., Walther, D. and 
#'  Selbig, J., 2007. \href{https://doi.org/10.1093/bioinformatics/btm069}{doi}.
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
#' pp <- ppcapM(t(missing.dataset), nPcs = 5)
#' names(pp)
#' 
#' # sigmasq estimation
#' abs(pp$sigmaSq-sigma.sq)
#' 
#' # X reconstruction
#' recon.X <- pp$pcaMethodsRes@loadings%*%t(pp$pcaMethodsRes@scores)
#' norm(recon.X-X, type="F")^2/(length(X))
#' 
#' # covariance estimation
#' norm(pp$Sigma-Sigma, type="F")^2/(length(X))
ppcapM <- function(myMat, nPcs=2, seed=NA, threshold=1e-4, maxIterations=1000,
                   loglike = TRUE, verbose=TRUE) {

  if (!is.na(seed)) 
    set.seed(seed)
  
  N <- nrow(myMat)
  D <- ncol(myMat)
  
  hidden   <- which(is.na(myMat))
  nMissing <- length(hidden)
  
  mu <- colMeans(myMat, na.rm = TRUE)

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
  Worth <- orthMat(ppcaOutput$W)
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
  pcaMethodsRes@missing  <- is.na(myMatsaved)
  
  # create hinton diagram
  if(verbose){
    plotrix::color2D.matplot(ppcaOutput$W,
                             extremes=c("black","white"),
                             main="Hinton diagram of loadings",
                             Hinton=TRUE)
  }
  
  if (loglike){
    # compute log-likelihood scores
    loglikeobs <- compute_loglikeobs(dat = t(myMatsaved), covmat = ppcaOutput$C,
                                     meanvec = mu,
                                     verbose = verbose)
    
    loglikeimp <- compute_loglikeimp(dat = t(myMatsaved), A = ppcaOutput$W,
                                     S = t(X),
                                     covmat = ppcaOutput$C, 
                                     meanvec = mu,
                                     verbose = verbose)
  }
  # Also return the standard ppcaNet output:
  output <- list()
  output[["W"]]              <- ppcaOutput$W
  output[["sigmaSq"]]        <- ppcaOutput$ss
  output[["Sigma"]]          <- ppcaOutput$C
  output[["m"]]             <- mu
  if (loglike) {
    output[["logLikeObs"]] <- loglikeobs
    output[["logLikeImp"]] <- loglikeimp
  }
  output[["pcaMethodsRes"]]  <- pcaMethodsRes
  
  return(output)
  
}


#' @title Calculate an orthonormal basis
#' 
#' @description A copied (unexported) function from \code{\link{pcaMethods}}. 
#'  ONB = orth(mat) is an orthonormal basis for the range of matrix mat. That is,
#'  ONB' * ONB = I, the columns of ONB span the same space as the columns of mat, 
#'  and the number of columns of ONB is the rank of mat.
#' 
#' @param mat \code{matrix} -- matrix to calculate the orthonormal basis of
#' @param skipInac \code{logical} -- do not include components with precision below
#'  \code{.Machine$double.eps} if \code{TRUE}
#'  
#' @return orthonormal basis for the range of \code{mat}
#' 
#' @seealso \code{\link[pcaMethods:orth]{orth}}
#' 
#' @author Wolfram Stacklies
#' 
#' @references Stacklies, W., Redestig, H., Scholz, M., Walther, D. and 
#'  Selbig, J., 2007. \href{https://doi.org/10.1093/bioinformatics/btm069}{doi}.
#'  
#' @examples 
#' set.seed(102)
#' X <- matrix(rnorm(10), 5, 2)
#' norm(X[,1], type="2")
#' norm(X[,2], type="2")
#' t(X[,1])%*%X[,2]
#' Xorth <- orthMat(X)
#' # now unit norms
#' norm(Xorth[,1], type="2")
#' norm(Xorth[,2], type="2")
#' # and zero dot product
#' t(Xorth[,1])%*%Xorth[,2]
orthMat <- function(mat, skipInac = FALSE) 
{
  if (nrow(mat) > ncol(mat)) {
    leftSVs <- ncol(mat)
  }
  else {
    leftSVs <- nrow(mat)
  }
  result <- svd(mat, nu = leftSVs, nv = ncol(mat))
  U <- result[[2]]
  S <- result[[1]]
  V <- result[[3]]
  m <- nrow(mat)
  n <- ncol(mat)
  if (m > 1) {
    s <- diag(S, nrow = length(S))
  }
  else if (m == 1) {
    s <- S[1]
  }
  else {
    s <- 0
  }
  tol <- max(m, n) * max(s) * .Machine$double.eps
  r <- sum(s > tol)
  if (r < ncol(U)) {
    if (skipInac) {
      warning("Precision for components ", r + 1, " - ", 
              ncol(U), " is below .Machine$double.eps. \n", 
              "Results for those components are likely to be inaccurate!!\n", 
              "These component(s) are not included in the returned solution!!\n")
    }
    else {
      warning("Precision for components ", r + 1, " - ", 
              ncol(U), " is below .Machine$double.eps. \n", 
              "Results for those components are likely to be inaccurate!!\n")
    }
  }
  if (skipInac) {
    ONB <- U[, 1:r, drop = FALSE]
    rownames(ONB) <- labels(mat[, 1:r, drop = FALSE])[[1]]
    colnames(ONB) <- labels(mat[, 1:r, drop = FALSE])[[2]]
  }
  else {
    ONB <- U
    rownames(ONB) <- labels(mat)[[1]]
    colnames(ONB) <- labels(mat)[[2]]
  }
  return(ONB)
}

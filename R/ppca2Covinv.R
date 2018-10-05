#' Inverse covariance matrix computation from PPCA
#'
#' @description Efficient inversion of the covariance matrix
#'   estimated from PPCA. 
#'   
#' @param ppcaOutput \code{list} -- the output object from running any
#'   of the PPCA functions in this package.
#'
#' @details The computation exploits the
#'   Woodbury identity so that a kxk matrix (where k is often less than 10)
#'   is inverted instead of the potentially large pxp matrix. The closed-form
#'   expression for the inverse depends upon parameters that are estimated
#'   in the PPCA algorithm.
#'
#' @return \code{matrix} -- the inverse of the covariance matrix.
#' @export
#'
#' @examples
ppca2Covinv <- function(ppcaOutput){
  # get the relevant parameters from the ppca output
  sigmaSq <- ppcaOutput$sigmaSq
  W <- ppcaOutput$W
  k <- ncol(W)
  p <- nrow(W)
  
  # compute the kxk matrix to be inverted and invert it
  M <- sigmaSq*diag(k)+crossprod(W)
  Minv <- chol2inv(chol(M))
  
  # compute the covariance inverse and return
  # the formula for the inversion comes from applying the
  # Woodbury formula to the covariance matrix decomposition
  # of ppca
  sigmaInv <-  (1/sigmaSq)*(diag(p)-W%*%tcrossprod(Minv,W))
  return(sigmaInv)
}
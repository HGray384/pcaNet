#' Internal cross-validation can be used for estimating the level of
#' structure in a data set and to optimise the choice of number of
#' principal components.
#'
#' A wrapper for the \code{\link[pcaMethods:Q2]{Q2}} function from
#' \code{\link{pcaMethods}}, which calculates \eqn{Q^2} for a PCA model.
#' This is the cross-validated version of \eqn{R^2} and can be interpreted as 
#' the ratio of variance that can be predicted independently by the PCA
#' model. Poor (low) \eqn{Q^2} indicates that the PCA model only
#' describes noise and that the model is unrelated to the true data
#' structure. The definition of \eqn{Q^2} is: \deqn{Q^2=1 -
#' \frac{\sum_{i}^{p}\sum_{j}^{n}(X -
#' \hat{X})^2}{\sum_{i}^{p}\sum_{j}^{n}X^2}}{Q^2=1 - sum_i^p
#' sum_j^n (X - \hat{X})^2 / \sum_i^p \sum_j^n(X^2)} for the matrix
#' \eqn{X} which has \eqn{n} rows and \eqn{p} columns. For a given
#' number of PC's X is estimated as \eqn{\hat{X}=TP'} (T are scores
#' and P are loadings). Although this defines the leave-one-out
#' cross-validation this is not what is performed if fold is less
#' than the number of rows and/or columns.  In 'impute' type CV,
#' diagonal rows of elements in the matrix are deleted and the
#' re-estimated.  In 'krzanowski' type CV, rows are sequentially left
#' out to build fold PCA models which give the loadings. Then,
#' columns are sequentially left out to build fold models for
#' scores. By combining scores and loadings from different models, we
#' can estimate completely left out values.  The two types may seem
#' similar but can give very different results, krzanowski typically
#' yields more stable and reliable result for estimating data
#' structure whereas impute is better for evaluating missing value
#' imputation performance. Note that since Krzanowski CV operates on
#' a reduced matrix, it is not possible estimate Q2 for all
#' components and the result vector may therefore be shorter than
#' \code{nPcs(object)}.
#' @title Cross-validation for PCA
#' @param obj A \code{pcaRes} object (result from previous PCA
#' analysis.)
#' @param originalData The matrix (or ExpressionSet) that used to
#' obtain the pcaRes object.
#' @param fold The number of groups to divide the data in.
#' @param nruncv The number of times to repeat the whole
#' cross-validation
#' @param type krzanowski or imputation type cross-validation
#' @param verbose \code{boolean} If TRUE Q2 outputs a primitive
#' progress bar.
#' @param variables indices of the variables to use during
#' cross-validation calculation. Other variables are kept as they are
#' and do not contribute to the total sum-of-squares.
#' @param ... Further arguments passed to the \code{\link{pca}} function called
#' within Q2.
#' @return A matrix or vector with \eqn{Q^2} estimates.
#' @export
#' @references Krzanowski, WJ. Cross-validation in principal
#' component analysis. Biometrics. 1987(43):3,575-584
#' @seealso \code{\link[pcaMethods:Q2]{Q2}}
#' @examples
#' # analogously to pcaMethods...
#' data(iris)
#' x <- iris[,1:4]
#' pcIr <- pcapM(as.matrix(x), nPcs=3, method="ppca", seed=104, scale="none")
#' q2 <- ppcaQ2(pcIr)
#' barplot(q2, main="Krzanowski CV", xlab="Number of PCs",
#'  ylab=expression(Q^2))
#' @author Henning Redestig, Ondrej Mikula
ppcaQ2 <- function (obj, 
                    originalData=obj$pcaMethodsRes@completeObs,
                    fold=5, nruncv=1, 
                    type=c("krzanowski", "impute"), verbose=interactive(),
                    variables=1:(obj$pcaMethodsRes@nVar), ...) 
{
 pcaMethods::Q2(object = obj$pcaMethodsRes, originalData=originalData,
                fold=fold, nruncv=nruncv, type=type, verbose=verbose,
                variables=variables, ...)
}

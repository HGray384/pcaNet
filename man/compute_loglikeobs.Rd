% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca_full.R
\name{compute_loglikeobs}
\alias{compute_loglikeobs}
\title{Compute the log-likelihood of the observed data given PCA parameter estimates}
\usage{
compute_loglikeobs(dat, covmat, meanvec, verbose = TRUE)
}
\arguments{
\item{dat}{\code{matrix} -- the data matrix with variables in rows and 
observations in columns.}

\item{covmat}{\code{matrix} -- the estimated covariance matrix.}

\item{meanvec}{\code{numeric} -- the estimated mean vector.}

\item{verbose}{\code{logical} -- whether extra output should be displayed.}
}
\value{
the log-likelihood value
}
\description{
The log-likelihood of the data for probabilistic PCA is known to be
 multivariate Gaussian. Using this, one can check the log-likelihood value of the
 observed data values given the parameter estimates from the PCA model. This can 
 be useful to compare different models.
}
\examples{
p <- 20
n <- 7
set.seed(10045)
X <- matrix(rnorm(p*n), p, n)
miss.inds <- sample(1:(p*n), (p*n)/4)
X[miss.inds] <- NA
M <- !is.na(X)
Nobs_i <- rowSums(M)
Mu <- rowSums(X, na.rm = TRUE) / Nobs_i
Mu2 <- rep(0, p)
covmat <- diag(p)
# using sample mean
compute_loglikeobs(dat=X, covmat=covmat, meanvec=Mu, verbose=TRUE)
# using zero mean
compute_loglikeobs(dat=X, covmat=covmat, meanvec=Mu2, verbose=TRUE)
}

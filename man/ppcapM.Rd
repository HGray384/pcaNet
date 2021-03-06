% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppcapM.R
\name{ppcapM}
\alias{ppcapM}
\title{Probabilistic PCA (pcaMethods version)}
\usage{
ppcapM(
  myMat,
  nPcs = 2,
  seed = NA,
  threshold = 1e-04,
  maxIterations = 1000,
  loglike = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{myMat}{\code{matrix} -- Pre-processed matrix (centered,
scaled) with variables in columns and observations in rows. The
data may contain missing values, denoted as \code{NA}.}

\item{nPcs}{\code{numeric} -- Number of components used for
re-estimation. Choosing few components may decrease the
estimation precision.}

\item{seed}{\code{numeric} -- the random number seed used, useful
to specify when comparing algorithms.}

\item{threshold}{\code{numeric} -- Convergence threshold. 
If the increase in precision of an update
falls below this then the algorithm is stopped.}

\item{maxIterations}{\code{numeric} -- Maximum number of estimation
steps.}

\item{loglike}{\code{logical} -- should the log-likelihood
of the estimated parameters be returned? See Details.}

\item{verbose}{\code{logical} -- verbose intermediary 
algorithm output.}
}
\value{
{A \code{list} of 4 elements:
\describe{
\item{W}{\code{matrix} -- the estimated loadings.}
\item{sigmaSq}{\code{numeric} -- the estimated isotropic variance.}
\item{Sigma}{\code{matrix} -- the estimated covariance matrix.}
\item{pcaMethodsRes}{\code{class} -- 
  see \code{\link[pcaMethods:pcaRes-class]{pcaRes}}.}
}}
}
\description{
Implements a probabilistic PCA missing value estimator, as in pcaMethods.
  Use of Rcpp makes this version faster and 
  the emphasised output is the covariance matrix \code{Sigma}, which 
  can be used for network reconstruction.
}
\details{
Details about the probabilistic model underlying PPCA are found in
Bishop 1999. The algorithm (Porta, 2005) uses an expectation maximisation
approach together with a probabilistic model to approximate the
principal axes (eigenvectors of the covariance matrix in PCA).
The estimation is done iteratively, the algorithm terminates if
either the maximum number of iterations is reached or if the
estimated increase in precision falls below \eqn{1e^{-4}}{1e^-4}.
}
\examples{
# simulate a dataset from a zero mean factor model X = Wz + epsilon
# start off by generating a random binary connectivity matrix
n.factors <- 5
n.genes <- 200
# with dense connectivity
# set.seed(20)
conn.mat <- matrix(rbinom(n = n.genes*n.factors,
                          size = 1, prob = 0.7), c(n.genes, n.factors))

# now generate a loadings matrix from this connectivity
loading.gen <- function(x){
  ifelse(x==0, 0, rnorm(1, 0, 1))
}

W <- apply(conn.mat, c(1, 2), loading.gen)

# generate factor matrix
n.samples <- 100
z <- replicate(n.samples, rnorm(n.factors, 0, 1))

# generate a noise matrix
sigma.sq <- 0.1
epsilon <- replicate(n.samples, rnorm(n.genes, 0, sqrt(sigma.sq)))

# by the ppca equations this gives us the data matrix
X <- W\%*\%z + epsilon
WWt <- tcrossprod(W)
Sigma <- WWt + diag(sigma.sq, n.genes)

# select 10\% of entries to make missing values
missFrac <- 0.1
inds <- sample(x = 1:length(X),
               size = ceiling(length(X)*missFrac),
               replace = FALSE)

# replace them with NAs in the dataset
missing.dataset <- X
missing.dataset[inds] <- NA

# run ppca
pp <- ppcapM(t(missing.dataset), nPcs = 5)
names(pp)

# sigmasq estimation
abs(pp$sigmaSq-sigma.sq)

# X reconstruction
recon.X <- pp$pcaMethodsRes@loadings\%*\%t(pp$pcaMethodsRes@scores)
norm(recon.X-X, type="F")^2/(length(X))

# covariance estimation
norm(pp$Sigma-Sigma, type="F")^2/(length(X))
}
\references{
Porta, J.M., Verbeek, J.J. and Kroese, B.J., 2005.
 \href{https://hal.inria.fr/inria-00321476/en}{link}

 Stacklies, W., Redestig, H., Scholz, M., Walther, D. and 
 Selbig, J., 2007. \href{https://doi.org/10.1093/bioinformatics/btm069}{doi}.
}

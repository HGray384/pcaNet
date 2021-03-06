% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppcapM.R
\name{orthMat}
\alias{orthMat}
\title{Calculate an orthonormal basis}
\usage{
orthMat(mat, skipInac = FALSE)
}
\arguments{
\item{mat}{\code{matrix} -- matrix to calculate the orthonormal basis of}

\item{skipInac}{\code{logical} -- do not include components with precision below
\code{.Machine$double.eps} if \code{TRUE}}
}
\value{
orthonormal basis for the range of \code{mat}
}
\description{
A copied (unexported) function from \code{\link{pcaMethods}}. 
 ONB = orth(mat) is an orthonormal basis for the range of matrix mat. That is,
 ONB' * ONB = I, the columns of ONB span the same space as the columns of mat, 
 and the number of columns of ONB is the rank of mat.
}
\examples{
set.seed(102)
X <- matrix(rnorm(10), 5, 2)
norm(X[,1], type="2")
norm(X[,2], type="2")
t(X[,1])\%*\%X[,2]
Xorth <- orthMat(X)
# now unit norms
norm(Xorth[,1], type="2")
norm(Xorth[,2], type="2")
# and zero dot product
t(Xorth[,1])\%*\%Xorth[,2]
}
\references{
Stacklies, W., Redestig, H., Scholz, M., Walther, D. and 
 Selbig, J., 2007. \href{https://doi.org/10.1093/bioinformatics/btm069}{doi}.
}
\seealso{
\code{\link[pcaMethods:orth]{orth}}
}
\author{
Wolfram Stacklies
}

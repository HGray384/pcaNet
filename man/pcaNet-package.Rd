\name{pcaNet-package}
\alias{pcaNet-package}
\alias{pcaNet}
\docType{package}
\title{
\packageTitle{pcaNet}
}
\description{
\packageDescription{pcaNet}
}
\details{
Algorithms for PPCA have been ported from the \code{PCAMV} MATLAB toolbox (Ilin and Raiko, 2010) and
extended from the \code{\link{pcaMethods}} (Stacklies et. al., 2007) R-package to focus on covariance matrix estimation and network reconstruction in the presence of missing values. Full PCA functionality with \code{\link{pcaMethods}} is retained in \code{pcaNet} due to the use of the \linkS4class{pcaRes} class.

The inverse of the covariance matrix from PPCA can be computed efficiently, and this functionality is provided in \code{\link{ppca2Covinv}}. Using the false discovery rate method from Strimmer (2008), the estimated partial correlations can be tested to construct a network. Whilst default behaviour for this is available, the full output of the testing is also provided, so that users may further explore the statistics using \code{\link{fdrtool}}. Functionality for visualising the covariance matrix is provided, as well as for the reconstructed network using \code{\link{igraph}} (Csardi and Nepusz, 2006).
}
\author{
\packageAuthor{pcaNet}

Maintainers: \packageMaintainer{pcaNet}
}
\references{
Oba, S., Sato, M.A., Takemasa, I., Monden, M.,
Matsubara, K.I. and Ishii, S., 2003. \href{https://doi.org/10.1093/bioinformatics/btg287}{doi}.

Stacklies, W., Redestig, H., Scholz, M., Walther, D. and 
Selbig, J., 2007. \href{https://doi.org/10.1093/bioinformatics/btm069}{doi}.

Ilin, A. and Raiko, T., 2010.
\href{http://www.jmlr.org/papers/v11/ilin10a.html}{link}

Porta, J.M., Verbeek, J.J. and Kroese, B.J., 2005.
\href{https://hal.inria.fr/inria-00321476/en}{link}

Strimmer, K., 2008.
\href{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-303}{link}.

Strimmer, K., 2008.
\href{https://doi.org/10.1093/bioinformatics/btn209}{doi}.

Csardi, G. and Nepusz, T., 2006.
\href{http://www.interjournal.org/manuscript_abstract.php?361100992.}{link}.
}
\keyword{ package }
\seealso{
\code{\link[pcaMethods:pcaMethods-package]{pcaMethods}}
}
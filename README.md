# pcaNet
Probabilistic principal components analysis for covariance matrix estimation and network reconstruction R-package (dev version)

Implementation of the algorithms from Gray and Kirk, (manuscript in progress) for probabilistic principal components 
analysis (PPCA) focussed on high-dimensional covariance matrix estimation and network reconstruction. PPCA
is a model-based approach for dimension reduction, which may be applied to incomplete data and implicitly performs 
covariance matrix estimation during model fitting. This package implements a number of PPCA algorithms and
provides extra functionality to visualise the estimated covariance matrix, efficiently compute its inverse, estimate
the underlying network of variables, and then also visualise this network.

See documentation files for information on the usage of PPCA.

This version can be directly installed from the R console using:

```
# install.packages("devtools")
devtools::install_github("HGray384/TAS")
```

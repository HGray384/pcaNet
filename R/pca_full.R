#' @title A wrapper for PCAMV (MATLAB) function implementations
#'
#' @description Implements the PPCA algorithms from (cite), previously
#'   only available in MATLAB. One element of the output is a pcaRes 
#'   object, providing an interface between PCAMV and pcaMethods.
#'
#' @param X \code{matrix} -- Data matrix with 
#'   observations in columns and variables in rows. The
#'   data may contain missing values, denoted as \code{NA},
#'   or \code{NaN}. 
#' @param ncomp \code{numeric} -- Number of components used for
#'   re-estimation. Choosing few components may decrease the
#'   estimation precision. Setting to \code{NA} results in 
#'   \code{ncomp} = min(n, p) -1, which will be slow for large
#'   data.
#' @param algorithm \code{c("ppca", "map", "vb")} -- the algorithm
#'   to be used for estimation, see Details. 
#' @param maxiters \code{numeric} -- Maximum number of estimation
#'   steps.
#' @param bias \code{logical} -- should the mean be estimated?
#' @param rotate2pca \code{logical} -- should the solution be rotated
#'   to a PCA basis? See Details.
#' @param loglike \code{logical} -- should the log-likelihood
#'   of the estimated parameters be returned? See Details.
#' @param verbose \code{logical} -- verbose intermediary 
#'   algorithm output.
#'
#' @details The \code{algorithm} argument provides the option of 
#'   performing either 'ppca' for PPCA, 'vb' for BPCA using a 
#'   variational approximation, or 'map' for a variational 
#'   approximation ignoring posterior uncertainty (for faster
#'   computation). See (cite) for the full models. Setting 
#'   \code{rotate2pca} will perform a post-estimation rotation of
#'   the scores and loadings matrices so that they satisfy the
#'   PCA conditions of orthonormality, see (cite) for the 
#'   derivations. \code{loglike} indicates whether 
#'   log-likelihood values for the resulting estimates should 
#'   be computed. This can be useful to compare different algorithms.
#'
#' @return {A \code{list} of 6 or 8 elements, depending on the value
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
#' \item{m}{\code{numeric} -- the number of iterations taken to 
#' converge.}
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
#' ppf <- pca_full(missing.dataset, ncomp=5, algorithm="vb", maxiters=5,
#' bias=TRUE, rotate2pca=FALSE, loglike=TRUE, verbose=TRUE)
pca_full <- function(X, ncomp=NA, algorithm = "vb", maxiters = 1000,
                     bias = TRUE, rotate2pca = TRUE, loglike = TRUE, 
                     verbose=TRUE){
  # comment this out before running
  # set.seed(20)
  # X <- missing.dataset
  # X <- matrix(rnorm(20000), 100, 200)
  # X[1, 2] <- NaN
  
  savedX <- X
  savedncomp <- ncomp
  savedalgorithm <- algorithm
  savedmaxiters <- maxiters
  savedbias <- bias
  savedrotate2pca = rotate2pca
  savedloglike = loglike
  savedverbose = verbose
  myCondNum <- 1000
  numOuterIts <- 1
  
  # while(myCondNum >= 1000 && numOuterIts < 6)
  # {
    X <- savedX
    ncomp <- savedncomp
    algorithm <- savedalgorithm
    maxiters <- savedmaxiters
    bias <- savedbias
    rotate2pca <- savedrotate2pca
    loglike <- savedloglike
    verbose <- savedverbose
    
    opts <- list(init='random',
                 maxiters=as.numeric(1000),
                 # bias=as.numeric(1),
                 # 'uniquesv'=0,
                 # 'autosave'=600,
                 # 'filename'='pca_f_autosave',
                 # 'minangle'=1e-8,
                 # 'algorithm'='vb',
                 niter_broadprior=as.numeric(100),
                 earlystop=as.numeric(0)
                 # 'rmsstop'=[ 100 1e-4 1e-3 ], # [] means no rms stop criteria
                 # 'cfstop'=[], # [] means no cost stop criteria
                 # 'verbose'=1,
                 # 'xprobe'=[],
                 # 'rotate2pca'=1,
                 # 'display'=0
    )
    
    if(algorithm == "vb")
    {
      use_prior = 1
      use_postvar = 1
    }
    else if(algorithm == "map")
    {
      use_prior = 1
      use_postvar = 0
      
    }
    else if(algorithm == "ppca")
    {
      use_prior = 0
      use_postvar = 0
    }
    else{
      stop("Algorithm must be one of 'vb', 'map' or 'ppca'")  
    }
    
    p <- nrow(X)
    n <- ncol(X)
    if (is.na(ncomp)){
      ncomp <- min(n, p)-1
    }
    
    # Missing values are marked as NaNs for consistency with MATLAB code,
    # from which this was translated
    if (any(is.na(X))){
      X[which(is.na(X))] <- NaN
    }
    M <- !is.nan(X)
    
    myMatsaved   <- X
    X[X==0]      <- .Machine$double.eps
    X[is.nan(X)] <- 0
    
    Nobs_i <- rowSums(M)
    
    
    notmiss <- which(X!=0, arr.ind = TRUE)
    IX      <- notmiss[,1]
    JX      <- notmiss[,2]
    ndata   <- length(IX)
    # clear data
    rm(notmiss)
    
    # % Compute indices Isv: Sv{Isv(j)} gives Sv for j, j=1...n2
    # these arguments are used to save memory when many Sv{j} are the same
    nobscomb <- n # not currently used but might be in a future release
    Isv      <- c()
    obscombj <- c()
    # end
    
    
    ####################################
    # parameter initialisation
    
    initialisedParms <- initParms(p, n, ncomp, verbose = verbose)
    A   <- initialisedParms$A
    S   <- initialisedParms$S
    Mu  <- initialisedParms$Mu
    V   <- initialisedParms$V
    Av  <- initialisedParms$Av
    Sv  <- initialisedParms$Sv
    Muv <- initialisedParms$Muv
    # 
    # write.table( A, file = "A.csv", row.names = F, col.names = F, sep = ",")
    # write.table( S, file = "S.csv", row.names = F, col.names = F, sep = ",")
    # write.table( Mu, file = "Mu.csv", row.names = F, col.names = F, sep = ",")
    # write.table( V, file = "V.csv", row.names = F, col.names = F, sep = ",")
    # write.table( Av, file = "Av.csv", row.names = F, col.names = F, sep = ",")
    # write.table( Sv, file = "Sv.csv", row.names = F, col.names = F, sep = ",")
    # write.table( Muv, file = "Muv.csv", row.names = F, col.names = F, sep = ",")
    
    
    ####################################
    
    Va  <- 1000*rep(1,ncomp)
    Vmu <- 1000
    
    if (!use_prior)
    {
      Va  <- Va *Inf
      Vmu <- Vmu *Inf 
    }
    
    # the use_postvar options from the matlab code are
    # executed in the pca_updates.cpp file linea 64-68
    if (!bias){
      # Muv is cleared in pca_updates.cpp line 72
      # Muv = c()
      Vmu = 0
    }
    
    if (is.null(Mu)){
      if (bias){
        Mu <- rowSums(X) / Nobs_i
      }else{
        Mu = rep(0, p)
      }
    }
    
    
    
    # data centering
    X <- subtractMu(Mu, X, M, p, n, bias, verbose = verbose) 
    ############################
    # compute initial rms
    ############################
    
    computedRMS <- compute_rms(X, A, S, M, ndata, verbose = verbose)
    errMx       <- computedRMS$errMx
    rms         <- computedRMS$rms
    
    ############################
    # initial hyperprior parameter values
    ############################
    
    hpVa <- 0.001
    hpVb <- 0.001
    hpV  <- 0.001
    
    #########################
    # CALL C++ FUNCTION
    #########################
    if (is.null(Isv)){Isv <- rep(0, 2)}
    IX <- IX -1 # C++ indexing
    JX <- JX -1 # C++ indexing
    if(verbose){
      verbose <- 1 # C++ true
    } else {
      verbose <- 0 # C++ false
    }
    if(bias){
      bias <- 1 # C++ true
    } else {
      bias <- 0 # C++ false
    }
    if(rotate2pca){
      rotate2pca <- 1 # C++ true
    } else {
      rotate2pca <- 0 # C++ false
    }
    
    # print('A init is \n')
    # print(A)
    # print('S init is \n')
    # print(S)
    # print('Mu init is \n')
    # print(Mu)
    # print('V init is \n')
    # print(V)
    # print('Av init is \n')
    # print(Av)
    # print('Sv init is \n')
    # print(Sv)
    # print('Muv init is \n')
    # print(Muv)
    # print('errMx init is \n')
    # print(errMx)
    # print('rms init is \n')
    # print(rms)
    
    ppcaOutput <- pca_updates(X=X, V=V, A=A, Va=Va, Av = Av, S = S, Sv = Sv, 
                              Mu = Mu, Muv = Muv, Vmu = Vmu,
                              hpVa = hpVa, hpVb = hpVb, hpV = hpV, ndata = ndata, Nobs_i = Nobs_i,
                              Isv = Isv, M = M, IX = IX, JX = JX, rms = rms, errMx = errMx, 
                              bias = bias, rotate2pca = rotate2pca, niter_broadprior = opts$niter_broadprior, 
                              use_prior = use_prior, use_postvar = use_postvar,
                              maxiters = maxiters, verbose = verbose)
    
    #########################
    # manage output
    #########################
    
    nPcs <- ncomp
    if(ppcaOutput$numIter == maxiters)
    {
      print('Maximum number of iterations reached')
    }
    
    bias <- as.logical(bias)
    
    # R2 computation
    R2cum      <- rep(NA, nPcs)
    TSS        <- sum(myMatsaved^2, na.rm = TRUE)
    
    for (i in 1:ncomp) {
      difference <- t(X) - (ppcaOutput$scores[, 1:i, drop = FALSE] %*% t(ppcaOutput$W[, 1:i, drop = FALSE]))
      R2cum[i] <- 1 - (sum(difference^2, na.rm = TRUE)/TSS)
    }
    
    
    pcaMethodsRes           <- new("pcaRes")
    pcaMethodsRes@scores    <- ppcaOutput$scores 
    pcaMethodsRes@loadings  <- ppcaOutput$W
    pcaMethodsRes@R2cum     <- R2cum
    pcaMethodsRes@method    <- algorithm
    pcaMethodsRes@missing   <- is.na(myMatsaved)
    pcaMethodsRes@nPcs <- ncomp # do we need to edit this?
    pcaMethodsRes@nObs <- ncol(myMatsaved)
    pcaMethodsRes@nVar <- nrow(myMatsaved)
    pcaMethodsRes@sDev <- apply(pcaMethodsRes@scores, 2, sd)
    pcaMethodsRes@center <- as.numeric(ppcaOutput$m)
    pcaMethodsRes@centered <- bias
    pcaMethodsRes@scale <- rep(1, nrow(myMatsaved))
    pcaMethodsRes@scaled <- "none"
    pcaMethodsRes@R2 <- pcaMethodsRes@R2cum[1]
    if (length(pcaMethodsRes@R2cum) > 1) {
      pcaMethodsRes@R2 <- c(pcaMethodsRes@R2,
                            diff(pcaMethodsRes@R2cum))
    }
    
    completeObs <- myMatsaved
    if(any(!M)){
      recData <- tcrossprod(pcaMethodsRes@scores[, 1:nPcs, drop = FALSE],
                            pcaMethodsRes@loadings[, 1:nPcs, drop = FALSE])
      # recData <- sweep(recData, 2, sc, "*")
      recData <- sweep(recData, 2, ppcaOutput$m, "+")
      completeObs[is.na(myMatsaved)] <- recData[is.na(myMatsaved)]
    }
    pcaMethodsRes@completeObs <- completeObs
    
    # create hinton diagram
    if(verbose){
      plotrix::color2D.matplot(ppcaOutput$W,
                               extremes=c("black","white"),
                               main="Hinton diagram of loadings",
                               Hinton=TRUE)
    }
    # back to R boolean
    if(verbose==1){
      verbose <- TRUE # C++ true
    } else {
      verbose <- FALSE # C++ false
    }
    
    if (loglike){
      # compute log-likelihood scores
      loglikeobs <- compute_loglikeobs(myMatsaved, ppcaOutput$C, ppcaOutput$m,
                                       verbose)
      
      loglikeimp <- compute_loglikeimp(myMatsaved, ppcaOutput$W, t(ppcaOutput$scores),
                                       ppcaOutput$C, ppcaOutput$m, verbose)
    }
    
    # Return standard ppcaNet output:
    
    output <- list()
    output[["W"]]              <- ppcaOutput$W
    output[["sigmaSq"]]        <- ppcaOutput$ss
    output[["Sigma"]]          <- ppcaOutput$C
    output[["m"]]              <- ppcaOutput$m
    if (loglike){
      output[["logLikeObs"]] <- loglikeobs
      output[["logLikeImp"]] <- loglikeimp
    }
    output[["numIter"]]  <- ppcaOutput$numIter
    output[["pcaMethodsRes"]]  <- pcaMethodsRes
    
    myCondNum   <- kappa(ppcaOutput$C, exact = TRUE)
    numOuterIts <- numOuterIts + 1
  # }
  # if (numOuterIts==6){
  #   warning("PPCA did not find a good solution after 5 runs.")
  # }
  return(output)
}



#' @title Initialise model parameters for \code{\link{pca_updates}}
#' 
#' @description Internal function within \code{\link{pca_full}} that initialises
#'  most model parameters. WARNING: does not initialise all parameters by itself
#'  correctly (since this depends on context) and so care should be taken when 
#'  using as a standalone function.
#'  
#' @details Random initialisations are set for the loadings and scores matrices.
#'  Diagonal matrices are set for the elements of \code{Av} and \code{Sv}. \code{V}
#'  is initialised to 1 and \code{Muv} is initialised to a vector of 1s.
#' 
#' @param p \code{numeric} -- the number of variables
#' @param n \code{numeric} -- the number of observations
#' @param ncomp \code{numeric} -- the number of components/latent variables
#' @param verbose \code{logical} -- whether extra output should be displayed
#' 
#' @return {A \code{list} of length 7:
#'  \describe{
#'  \item{A}{\code{matrix} -- initialised loadings matrix with observed variables
#'   in rows and latent variables in columns.}
#'  \item{S}{\code{matrix} -- initialised factor scores matrix with latent variables
#'   in rows and observations in columns.}
#'  \item{V}{\code{numeric} -- scalar value corresponding to the initialised 
#'  variance of the error parameter.}
#'  \item{Av}{\code{array} -- initialised covariance matrices of the rows of A.}
#'  \item{Sv}{\code{array} -- initialised covariance matrices of the rows of S.}
#'  \item{Muv}{\code{numeric} -- the initialisation of the prior variance of Mu.}
#'  }
#' }
#' 
#' @examples 
#' init.model <- initParms(p=10, n=10, ncomp=2, verbose = TRUE)
#' init.model$A
#' init.model$Av
initParms <- function(p, n, ncomp, verbose=TRUE)
{
  if(verbose) {
    cat("Initialising parameters... \n")
  }
  randNumMatrix <- matrix(rnorm(p*ncomp), nrow = p, ncol = ncomp)
  qrDecomp      <- qr(randNumMatrix)
  # if ncomp > p then need to take R instead of Q
  if (ncomp < p) {
    A             <- qr.Q(qrDecomp)
  } else {
    A             <- qr.R(qrDecomp)
  }
  # A <- matrix(1, nrow = p, ncol = ncomp)
  
  Av <- lapply(1:p, function(x){diag(ncomp)})
  Av <- simplify2array(Av)
  
  # Mu = m
  Mu <- c()
  Muv <- rep(1, p)
  
  # V = vy
  V <- 1
  
  # S = X
  S  <- matrix(rnorm(ncomp*n), nrow = ncomp, ncol = n)
  # S  <- matrix(1, nrow = ncomp, ncol = n)
  Sv <- lapply(1:n, function(x){diag(ncomp)})
  Sv <- simplify2array(Sv)
  return(list(A = A, S = S, Mu = Mu, V = V, Av = Av, Sv = Sv, Muv = Muv))
}

#' @title Subtract the row means from a matrix of data with missing values
#' 
#' @description internal function within \code{\link{pca_full}} to subtract the row means
#'  from a matrix of data using only the observed values. Offers little utility 
#'  standalone.
#' 
#' @param Mu \code{numeric} -- the sample mean of the observed variables.
#' @param X \code{matrix} -- the data matrix with variables in rows and 
#'  observations in columns.
#' @param M \code{matrix} -- logical matrix whose values indicate whether
#'  the corresponding entry in \code{X} is observed.
#' @param p \code{numeric} -- the number of variables.
#' @param n \code{numeric} -- the number of observations.
#' @param update_bias \code{logical} -- whether the mean should be subtracted.
#'  or not.
#' @param verbose \code{logical} -- whether extra output should be displayed.
#' 
#' @return \code{X} \code{matrix} -- centered data matrix.
#' 
#' @examples
#' p <- 20
#' n <- 7
#' set.seed(10045)
#' X <- matrix(rnorm(p*n), p, n)
#' miss.inds <- sample(1:(p*n), (p*n)/4)
#' X[miss.inds] <- NA
#' M <- !is.na(X)
#' Nobs_i <- rowSums(M)
#' Mu <- rowSums(X, na.rm = TRUE) / Nobs_i
#' update_bias <- TRUE
#' Xcent <- subtractMu(Mu=Mu, X=X, M=M, p=p, n=n, update_bias=update_bias, verbose=TRUE)
#' X-Xcent
#' Mu # all observed values in each column equal to Mu
subtractMu <- function(Mu, X, M, p, n, update_bias, verbose=TRUE)
{
  if(verbose) {
    cat("Centering data... \n")
  }
  if (update_bias){
    X <- X - matrix(Mu,p,n)*M; 
  }
  return(X)
}


#' @title Compute the root mean-squared error of a PCA projection
#' 
#' @description Root mean-squared error is the square root of the element-wise error's mean.
#'  This is a useful quantity to display during parameter estimation in \code{\link{pca_updates}}
#'  since it is a measure of how well the PCA projection is fitting the data.
#' 
#' @param X \code{matrix} -- the data matrix with variables in rows and 
#'  observations in columns.
#' @param A \code{matrix} -- initialised loadings matrix with observed variables
#'   in rows and latent variables in columns.
#' @param S \code{matrix} -- initialised factor scores matrix with latent variables
#'   in rows and observations in columns.
#' @param M \code{matrix} -- logical matrix whose values indicate whether
#'  the corresponding entry in \code{X} is observed.
#' @param ndata \code{numerical} -- the total number of observed values.
#' @param verbose \code{logical} -- whether extra output should be displayed.
#' 
#' @return {A \code{list} of length 2:
#'  \describe{
#'  \item{errMx}{\code{matrix} -- matrix of element-wise differences (errors)
#'  between the observed data and the PCA projection.}
#'  \item{rms}{\code{numerical} -- root mean-squared error of the PCA
#'  projection.}
#'  }
#' }
#' 
#' @examples 
#' p <- 20
#' n <- 7
#' set.seed(10045)
#' X <- matrix(rnorm(p*n), p, n)
#' miss.inds <- sample(1:(p*n), (p*n)/4)
#' X[miss.inds] <- NA
#' M <- !is.na(X)
#' Nobs_i <- rowSums(M)
#' Mu <- rowSums(X, na.rm = TRUE) / Nobs_i
#' update_bias <- TRUE
#' Xcent <- subtractMu(Mu=Mu, X=X, M=M, p=p, n=n, update_bias=update_bias, verbose=TRUE)
#' init.model <- initParms(p=p, n=n, ncomp=2, verbose = TRUE)
#' compute_rms(X=X, A=init.model$A, S=init.model$S, M=M, ndata=sum(Nobs_i), verbose=TRUE)
compute_rms <- function(X, A, S, M, ndata, verbose=TRUE)
{
  if(verbose) {
    cat("Computing rms... \n")
  }
  errMx <- (X - A%*%S)*M
  rms   <- sqrt(sum(errMx^2, na.rm = TRUE)/ndata)
  list(errMx = errMx, rms = rms)
}


#' @title Compute the log-likelihood of the observed data given PCA parameter estimates
#' 
#' @description The log-likelihood of the data for probabilistic PCA is known to be
#'  multivariate Gaussian. Using this, one can check the log-likelihood value of the
#'  observed data values given the parameter estimates from the PCA model. This can 
#'  be useful to compare different models.
#' 
#' @param dat \code{matrix} -- the data matrix with variables in rows and 
#'  observations in columns.
#' @param covmat \code{matrix} -- the estimated covariance matrix.
#' @param meanvec \code{numeric} -- the estimated mean vector.
#' @param verbose \code{logical} -- whether extra output should be displayed.
#' 
#' @return the loglikelihood value
#' 
#' @examples 
#' p <- 20
#' n <- 7
#' set.seed(10045)
#' X <- matrix(rnorm(p*n), p, n)
#' miss.inds <- sample(1:(p*n), (p*n)/4)
#' X[miss.inds] <- NA
#' M <- !is.na(X)
#' Nobs_i <- rowSums(M)
#' Mu <- rowSums(X, na.rm = TRUE) / Nobs_i
#' Mu2 <- rep(0, p)
#' covmat <- diag(p)
#' # using sample mean
#' compute_loglikeobs(dat=X, covmat=covmat, meanvec=Mu, verbose=TRUE)
#' # using zero mean
#' compute_loglikeobs(dat=X, covmat=covmat, meanvec=Mu2, verbose=TRUE)
compute_loglikeobs <- function(dat, covmat, meanvec, verbose=TRUE)
{
  if (verbose){
    cat("Computing log-likelihood of observed data... \n")
  }
  # convert any NAs to NaNs for convenience
  tmp <- dat
  if (any(is.na(dat))){
    tmp[which(is.na(dat))] <- NaN
  }
  # check if there are missing values
  if(any(is.nan(tmp))){
    # check which samples contain missing values
    missing.cols <- which(is.nan(tmp), arr.ind = TRUE)[,2]
    if(length(missing.cols)>1){
      # if there are multiple samples with missing values
      # then we can 'apply' over them
      loglikemiss <- apply(tmp[,missing.cols], 2, function(xcol.elems){
        missing.inds <- which(is.nan(xcol.elems));
        x.obs.tmp <- xcol.elems[-missing.inds];
        mu.obs.tmp <- meanvec[-missing.inds];
        cov.obs.tmp <- covmat[-missing.inds, -missing.inds];
        dens <- mvtnorm::dmvnorm(x = x.obs.tmp, mean = mu.obs.tmp, 
                                 sigma = cov.obs.tmp, log = TRUE);
        dens
      })
    } else {
      # if there is only one sample with missing values
      # then we have no need for an 'apply'
      xcol.elems <- tmp[,missing.cols]
      missing.inds <- which(is.nan(xcol.elems));
      x.obs.tmp <- xcol.elems[-missing.inds];
      mu.obs.tmp <- meanvec[-missing.inds];
      cov.obs.tmp <- covmat[-missing.inds, -missing.inds];
      dens <- mvtnorm::dmvnorm(x = x.obs.tmp, mean = mu.obs.tmp, 
                               sigma = cov.obs.tmp, log = TRUE);
      loglikemiss <- dens
    }
    # loglikelihood for samples with no missing values
    if ((length(missing.cols)+1)<ncol(tmp)) {
      # multiple samples with no missing values
      loglikenotmiss <- apply(tmp[,-missing.cols], 2, function(xcol.full){
        dens2 <- mvtnorm::dmvnorm(x = xcol.full, mean = meanvec, 
                                  sigma = covmat, log = TRUE);
        dens2
      })
    } else if (length(missing.cols)==(ncol(tmp)-1)){
      # a single sample with no missing values
      loglikenotmiss <- mvtnorm::dmvnorm(x = tmp[,-missing.cols], 
                                         mean = meanvec, 
                                         sigma = covmat, 
                                         log = TRUE);
    } else {
      # all samples have missing values
      loglikenotmiss <- c()
    }
    loglike <- c(loglikemiss, loglikenotmiss)
  } else {
    # no missing values
    loglike <- apply(tmp, 2, function(xcol.elems){
      dens <- mvtnorm::dmvnorm(x = xcol.elems, mean = meanvec, 
                               sigma = covmat, log = TRUE);
      dens
    })
  }
  return(sum(loglike))
}


#' @title Compute the log-likelihood of the observed data given PCA parameter estimates
#' 
#' @description The log-likelihood of the data for probabilistic PCA is known to be
#'  multivariate Gaussian. Using this, one can check the log-likelihood value of the
#'  observed data values given the parameter estimates from the PCA model. This can 
#'  be useful to compare different models.
#' 
#' @param dat \code{matrix} -- the data matrix with variables in rows and 
#'  observations in columns.
#' @param A \code{matrix} -- estimated loadings matrix with observed variables
#'   in rows and latent variables in columns.
#' @param S \code{matrix} -- estimated factor scores matrix with latent variables
#'   in rows and observations in columns.
#' @param covmat \code{matrix} -- the estimated covariance matrix.
#' @param meanvec \code{numeric} -- the estimated mean vector.
#' @param verbose \code{logical} -- whether extra output should be displayed.
#' 
#' @return the loglikelihood value
#' 
#' @examples 
#' p <- 20
#' n <- 20
#' set.seed(10045)
#'   verbose <- 1
#'   bias <- 1
#'   rotate2pca <- 1
#'   ncomp <- 2
#'   maxiters <- 1000
#'   opts <- list(init='random',
#'    maxiters=as.numeric(1000),
#'    niter_broadprior=as.numeric(100),
#'    earlystop=as.numeric(0)
#'   )
#'   use_prior = 1
#'   use_postvar = 1
#'   X <- matrix(rnorm(p*n), p, n)
#'   miss.inds <- sample(1:(p*n), round(p*n/10))
#'   X[miss.inds] <- NaN
#'   Xsaved <- X
#'   M <- !is.nan(X)
#'   X[X==0]      <- .Machine$double.eps
#'   X[is.nan(X)] <- 0
#'   
#'   notmiss <- which(X!=0, arr.ind = TRUE)
#'   IX      <- notmiss[,1]
#'   JX      <- notmiss[,2]
#'   
#'   Nobs_i = rowSums(M)
#'   ndata   <- length(IX)
#'   # C++ indexing
#'   IX <- IX -1 
#'   JX <- JX -1
#'   
#'   initialisedParms <- initParms(p, n, ncomp, verbose = verbose)
#'   A   <- initialisedParms$A
#'   S   <- initialisedParms$S
#'   Mu  <- initialisedParms$Mu
#'   V   <- initialisedParms$V
#'   Av  <- initialisedParms$Av
#'   Sv  <- initialisedParms$Sv
#'   Muv <- initialisedParms$Muv
#'   Va  <- 1000*rep(1,ncomp)
#'   Vmu <- 1000
#'   Mu <- rowSums(X) / Nobs_i
#'   computedRMS <- compute_rms(X, A, S, M, ndata, verbose = verbose)
#'   errMx       <- computedRMS$errMx
#'   rms         <- computedRMS$rms
#'   hpVa <- 0.001
#'   hpVb <- 0.001
#'   hpV  <- 0.001
#'   Isv <- rep(0, 2)
#'   # data centering
#'   X <- subtractMu(Mu, X, M, p, n, bias, verbose = verbose) 
#'   ppcaOutput <- pca_updates(X=X, V=V, A=A, Va=Va, Av = Av, S = S, Sv = Sv, 
#'   Mu = Mu, Muv = Muv, Vmu = Vmu,
#'   hpVa = hpVa, hpVb = hpVb, hpV = hpV, ndata = ndata, Nobs_i = Nobs_i,
#'   Isv = Isv, M = M, IX = IX, JX = JX, rms = rms, errMx = errMx, 
#'   bias = bias, rotate2pca = rotate2pca, niter_broadprior = opts$niter_broadprior, 
#'   use_prior = use_prior, use_postvar = use_postvar,
#'   maxiters = maxiters, verbose = verbose)
#'   # initialised model log-likelihood
#'   compute_loglikeimp(dat=Xsaved, A=A, S=S, covmat=tcrossprod(A)+diag(p),
#'   meanvec=Mu, verbose=TRUE)
#'   # estimated model log-likelihood
#'   compute_loglikeimp(dat=Xsaved, A=ppcaOutput$W, S=t(ppcaOutput$scores), covmat=ppcaOutput$C,
#'   meanvec=ppcaOutput$m, verbose=TRUE)
compute_loglikeimp <- function(dat, A, S, covmat, meanvec, verbose=TRUE)
{
  if (verbose) {
    cat("Computing log-likelihood of data with imputed missing values... \n")
  }
  # convert any NAs to NaNs for convenience
  tmp <- dat
  if (any(is.na(dat))){
    tmp[which(is.na(dat))] <- NaN
  }
  if (!any(is.nan(tmp))){
    return("There were no missing values to impute")
  }
  proj <- A%*%S
  tmp[which(is.nan(tmp))] <- proj[which(is.nan(tmp))]
  loglike <- apply(tmp, 2, function(xcol.elems){
    x.imp.tmp <- xcol.elems;
    dens <- mvtnorm::dmvnorm(x = x.imp.tmp, mean = meanvec, 
                             sigma = covmat, log = TRUE);
    dens
  })
  rm(tmp)
  rm(proj)
  return(sum(loglike))
}
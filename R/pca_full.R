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
  
  while(myCondNum >= 1000 && numOuterIts < 5)
  {
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
    M = !is.nan(X)
    
    myMatsaved   <- X
    X[X==0]      <- .Machine$double.eps
    X[is.nan(X)] <- 0
    
    Nobs_i = rowSums(M)
    
    
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
  }
  if (numOuterIts==6){
    warning("PPCA did not find a good solution after 5 runs.")
  }
  return(output)
}



# Initialise model parameters:
# ============================
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
  
  Av <- lapply(1:p, function(x){diag(ncomp)})
  Av <- simplify2array(Av)
  
  # Mu = m
  Mu <- c()
  Muv <- rep(1, p)
  
  # V = vy
  V <- 1
  
  # S = X
  S  <- matrix(rnorm(ncomp*n), nrow = ncomp, ncol = n)
  Sv <- lapply(1:n, function(x){diag(ncomp)})
  Sv <- simplify2array(Sv)
  return(list(A = A, S = S, Mu = Mu, V = V, Av = Av, Sv = Sv, Muv = Muv))
}


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



compute_rms <- function(X, A, S, M, ndata, verbose=TRUE)
{
  if(verbose) {
    cat("Computing rms... \n")
  }
  errMx <- (X - A%*%S)*M
  rms   <- sqrt(sum(errMx^2)/ndata)
  list(errMx = errMx, rms = rms)
}



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
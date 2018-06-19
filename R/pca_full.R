pca_full <- function(X, ncomp=NA, algorithm = "vb", maxiters = 1000, verbose=TRUE){
# comment this out before running
# set.seed(20)
# X <- missing.dataset
# X <- matrix(rnorm(20000), 100, 200)
# X[1, 2] <- NaN

opts <- list(init='random',
         maxiters=as.numeric(1000),
         bias=as.numeric(1),
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
  ncomp <- n-1
}

# Missing values are marked as NaNs
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

initialisedParms <- initParms(p, n, ncomp)
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

if (!opts$bias){
  Muv = c()
  Vmu = 0
}

if (is.null(Mu)){
  if (opts$bias){
    Mu <- rowSums(X) / Nobs_i
  }else{
    Mu = rep(0, p)
  }
}



# data centering
X <- subtractMu(Mu, X, M, p, n, opts$bias) 
############################
# compute initial rms
############################

computedRMS <- compute_rms(X, A, S, M, ndata)
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
ppcaOutput <- pca_updates(X=X, V=V, A=A, Va=Va, Av = Av, S = S, Sv = Sv, 
                          Mu = Mu, Muv = Muv, Vmu = Vmu,
         hpVa = hpVa, hpVb = hpVb, hpV = hpV, ndata = ndata, Nobs_i = Nobs_i,
         Isv = Isv, M = M, IX = IX, JX = JX, rms = rms, errMx = errMx, 
         bias = opts$bias, niter_broadprior = opts$niter_broadprior, 
         use_prior = use_prior, use_postvar = use_postvar,
         maxiters = maxiters, verbose = verbose)

#########################
# manage output
#########################

nPcs <- ncomp
if(ppcaOutput$numIter == opts$maxiters)
{
  print('Maximum number of iterations reached')
}
R2cum      <- rep(NA, nPcs)
TSS        <- sum(myMatsaved^2, na.rm = TRUE)

# usually compute the above variables but currently
# deciding whether to exclude shrunk dimensions

pcaMethodsRes           <- new("pcaRes")
pcaMethodsRes@scores    <- ppcaOutput$scores 
pcaMethodsRes@loadings  <- ppcaOutput$W
pcaMethodsRes@R2cum     <- rep(1, ncomp)
pcaMethodsRes@method    <- algorithm

# create hinton diagram
if(verbose){
plotrix::color2D.matplot(ppcaOutput$W,
                                  extremes=c("black","white"),
                                  main="Hinton diagram (white +, black -)",
                                  Hinton=TRUE)
}

# Return standard ppcaNet output:

output <- list()
output[["W"]]              <- ppcaOutput$W
output[["sigmaSq"]]        <- ppcaOutput$ss
output[["Sigma"]]          <- ppcaOutput$C
output[["numIter"]]  <- ppcaOutput$numIter
output[["pcaMethodsRes"]]  <- pcaMethodsRes

return(output)
}



# Initialise model parameters:
# ============================
initParms <- function(p, n, ncomp)
{
  randNumMatrix <- matrix(rnorm(p*ncomp), nrow = p, ncol = ncomp) 
  qrDecomp      <- qr(randNumMatrix)
  A             <- qr.Q(qrDecomp)
  
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


subtractMu <- function(Mu, X, M, p, n, update_bias)
{
  if (update_bias){
    X <- X - matrix(Mu,p,n)*M; 
  }
  return(X)
}



compute_rms <- function(X, A, S, M, ndata)
{
  errMx <- (X - A%*%S)*M
  rms   <- sqrt(sum(errMx^2)/ndata)
  list(errMx = errMx, rms = rms)
}
# %  PCA_FULL - PCA with unrestricted Gaussian posterior
# %
# %  PCA with unrestricted Gaussian pdfs (full covariance matrices) for
# %  approximating the posterior distributions of A(i,:) and S(:,j) in
# %  the model X(:,j) = Mu + A*S(:,j) + Noise. The noise is isotropic
# %  and Gaussian with the variance V.
# %
# %  [ A, S, Mu, V, CV, HP, LC ] = PCA_FULL( X, N ) identifies the model
# %  for the given data matrix X and number of principal components N.
# %  Matrix X can be either sparse with only observed values included
# %  (note that observed zeros should be replaced with eps) or a full
# %  matrix with missing values replaced by NaNs.
# %
# %  The function returns the mean values of the model parameters in A,
# %  S, Mu, and V. The posterior variances are stored in CV such that
# %  CV.A{j} is the posterior covariance matrix for A(i,:) and
# %  CV.S{CV.Isv(j)} (or CV.S{j} if CV.Isv is empty) is the same for
# %  S(:,j). HP contains point estimates of the hyperparameters: HP.Va
# %  is the prior variance for A(:,k) and HP.Vmu is the same for Mu.
# %  
# %  LC is a structure with learning curves (rms training and probing
#                                            %  errors, cost function values, time).
# %
# %  PCA_FULL( X, N, 'algorithm', name ) specifies the algorithm:
#   %   'ppca': Probabilistic PCA (no prior and point estimates for A, Mu)
# %   'map':  Prior on A(:,k) and Mu, point estimates for A(i,:) and Mu
# %   'vb':   Variational Bayesian PCA (prior on A and Mu, Gaussian
#                                       %           posterior approximation for A(i,:) and Mu) {default}
# %  
# %  Other optional parameter/value pairs with {default values}:
#   %   init       - Initialization type ({'random'} - random,
#                                         %                filename: load from file, structure: from given data)
# %   maxiters   - Maximum number of iterations {1000}
# %   rotate2pca - Whether to perform rotation of A and S to speed-up
# %                convergence {1}
# %   uniquesv   - Whether to compute only unique covariance matrices of
# %                S(:,j) {1}
# %   minangle   - Termination by minimum angle between subspaces
# %                defined by A {1e-8}
# %   rmsstop    - Termination by rms training error ([] - no rms stop)
# %                {[ 100 1e-4 1e-3 ]}: 100 iterations, 1e-4 absolute
# %                tolerance, 1e-3 relative tolerance,
# %   cfstop     - Termination by cost function value {[]} (similarly to
#                                                           %                rmsstop). The cost function is computed only if this
# %                option is nonempty
# %   xprobe     - Validation data set (of the same dimensions as X)
# %   earlystop  - Whether to use early stopping based on probing error
# %   verbose    - Progress information display level (0,{1},2)
# %   display    - Plot progress {0}
# %   autosave   - Auto-save after each {600} seconds
# %   filename   - Name of the file for auto-save {'pca_f_autosave'}
# %
# %  OUT = PCA_FULL( X, N ) returns all the outputs in a single
# %  structure variable OUT. Learning can be continued as follows:
#   %    out = pca_full( X, N );
#   %    [ A, S, Mu, V, CV, HP, LC ] = pca_full( X, N, 'init', out );
#   
#   %  TODO: More efficient implementation of V
#   
#   %  This software is provided "as is", without warranty of any kind.
#   %  Alexander Ilin, Tapani Raiko
  
pca_full <- function(X, ncomp=NA){
# comment this out before running
# set.seed(20)
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

# [ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
# if ~isempty(errmsg), error( errmsg ), end
# if ~isempty(wrnmsg), warning( wrnmsg ), end
# Xprobe = opts.xprobe;

# switch opts.algorithm
# case 'ppca'
# use_prior = 0;
# use_postvar = 0;
# case 'map'
# use_prior = 1;
# use_postvar = 0;
# case 'vb'
# use_prior = 1;
# use_postvar = 1;
# otherwise
# error( 'Wrong value of the argument ''algorithm''' )
# end

# only use VB options
# use_prior <- 1
# use_postvar <- 1

# 
# [n1x,n2x] = size(X);
# [ X, Xprobe, Ir, Ic, opts.init ] = rmempty( X, Xprobe,...
#                                             opts.init, opts.verbose );

p <- nrow(X)
n <- ncol(X)
if (is.na(ncomp)){
  ncomp <- n-1
}

# if issparse(X)
# % X is a sparse matrix with only observed values
# M = spones(X);    Mprobe = spones(Xprobe);
# else
#   % Missing values are marked as NaNs
M = !is.nan(X)
# Mprobe = ~isnan(Xprobe)
myMatsaved <- X
X[X==0] = .Machine$double.eps
# Xprobe(Xprobe==0) = eps;
X[is.nan(X)] = 0
# Xprobe(isnan(Xprobe)) = 0;
# 
# end

Nobs_i = rowSums(M)

# nprobe = nnz(Mprobe);
# if nprobe == 0
# Xprobe = [];
# opts.earlystop = 0;
# end
notmiss <- which(X!=0, arr.ind = TRUE)
IX <- notmiss[,1]
JX <- notmiss[,2]
ndata <- length(IX)
# missingI <- which(1:nrow(X)!=IX)
# missingJ <- which(1:ncol(X)!=JX)
# # [IX,JX,data] = find(X);
# clear data
rm(notmiss)

# % Compute indices Isv: Sv{Isv(j)} gives Sv for j, j=1...n2
# if opts.uniquesv
# [ nobscomb, obscombj, Isv ] = miscomb(M,opts.verbose);
# else
nobscomb <- n
Isv <- c()
obscombj <- c()
# end

# [ A, S, Mu, V, Av, Sv, Muv ] = ...
# InitParms( opts.init, n1, n2, ncomp, nobscomb, Isv );

####################################
####################################
# parameter initialisation
# A = orth(randn(p,ncomp));
# A = W
A <- qr.Q(qr(matrix(rnorm(p*ncomp), nrow = p, ncol = ncomp)))
Av <- lapply(1:p, function(x){diag(ncomp)})
Av <- simplify2array(Av)

# Mu = m
Mu <- c()
Muv <- rep(1, p)

# V = vy
V <- 1

# S = X
S <- matrix(rnorm(ncomp*n), nrow = ncomp, ncol = n)
Sv <- lapply(1:n, function(x){diag(ncomp)})
Sv <- simplify2array(Sv)

# Sv = sigmax
# if isfield( init, 'Sv' ) && ~isempty(init.Sv)
#   if nobscomb < n2
#     [B,I] = unique(Isv,'first');
#   if ~iscell(init.Sv)
#     Sv = cell(1,nobscomb);
#     for j = 1:nobscomb
#       Sv{j} = diag(init.Sv(:,Isv(I(j))));
# end
# elseif isfield( init, 'Isv' ) && ~isempty(init.Isv)
#   Sv = { init.Sv{ init.Isv(I) } };
# else
#   for j = 1:nobscomb
#     Sv{j} = init.Sv{Isv(I(j))};
#   end
# end
# else
#   if ~iscell(init.Sv)
#   Sv = cell(1,n2);
# for j = 1:n2, Sv{j} = diag(init.Sv(:,j)); end
#   
# elseif isfield( init, 'Isv' ) && ~isempty(init.Isv)
#   Sv = { init.Sv{ init.Isv } };
# elseif length(init.Sv) == n2
#   Sv = init.Sv;
# end
# end
# else
#   Sv = cell(1,nobscomb);
# for j = 1:nobscomb, Sv{j} = eye(ncomp); end
# end
####################################
####################################

# if use_prior
Va = 1000*rep(1,ncomp)
Vmu = 1000
# else
#   Va = repmat(inf,1,ncomp); Vmu = inf;
# end

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

# [X,Xprobe] = SubtractMu( Mu, X, M, Xprobe, Mprobe, opts.bias );
############################
############################

if (opts$bias){
  X <- X - matrix(Mu,p,n)*M; 
}

############################
############################

# [rms,errMx] = compute_rms( X, A, S, M, ndata );
# prms = compute_rms( Xprobe, A, S, Mprobe, nprobe );
errMx <- (X - A%*%S)*M
rms <- sqrt(sum(errMx^2)/ndata)

# lc.rms = rms; lc.prms = prms; lc.time = 0; lc.cost = NaN;

# dsph = DisplayInit( opts.display, lc );
# PrintFirstStep( opts.verbose, rms, prms );
# Aold <- A

# % Parameters of the prior for variance parameters
hpVa = 0.001
hpVb = 0.001
hpV = 0.001

# time_start = clock;
# time_autosave = time_start;
# tic
# %opts.niter_broadprior = 100;

#########################
# CALL C++ FUNCTION
#########################
if (is.null(Isv)){Isv <- rep(0, 2)}
IX <- IX -1 # C++ indexing
JX <- JX -1 # C++ indexing
ppcaOutput <- pca_updates(X=X, V=V, A=A, Va=Va, Av = Av, S = S, Sv = Sv, Mu = Mu, Muv = Muv, Vmu = Vmu,
         hpVa = hpVa, hpVb = hpVb, hpV = hpV, ndata = ndata, Nobs_i = Nobs_i, Isv = Isv,
         M = M, IX = IX, JX = JX, rms = rms, errMx = errMx, bias = opts$bias,
         niter_broadprior = opts$niter_broadprior, use_prior = 1, maxiters = 1000)

#########################


  
  
#   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   % Find the PCA rotation: This has to be checked
#   function [ dMu, A, Av, S, Sv ] = ...
#   RotateToPCA( A, Av, S, Sv, Isv, obscombj, update_bias );
#   
#   n1 = size(A,1);
#   n2 = size(S,2);
#   
#   if update_bias
#   mS = mean(S,2);
#   dMu = A*mS;
#   S = S - repmat(mS,1,n2);
#   else
#   dMu = 0;
#   end
#   
#   covS = S*S';
# if isempty(Isv)
# for j = 1:n2
# covS = covS + Sv{j};
# end
# else
#   nobscomb = length(obscombj);
# for j = 1:nobscomb
# covS = covS + ( length(obscombj{j})*Sv{j} );
# end
# end
# 
# covS = covS / n2;
# %covS = covS / (n2-n1);
# [VS,D] = eig(covS);
# RA = VS*sqrt(D);
# A = A*RA;
# covA = A'*A;
#   if ~isempty(Av)
#   for i = 1:n1
#   Av{i} = RA'*Av{i}*RA;
# covA = covA + Av{i};
# end
# end
# covA = covA / n1;
# [VA,DA] = eig(covA);
# [DA,I] = sort( -diag(DA) );
# DA = -DA;
# VA = VA(:,I);
# A = A*VA;
# 
# if ~isempty(Av)
# for i = 1:n1
# Av{i} = VA'*Av{i}*VA;
#   end
#   end
#   R = VA'*diag(1./sqrt(diag(D)))*VS';
#   
#   S = R*S;
#   for j = 1:length(Sv)
#   Sv{j} = R*Sv{j}*R';
# end
nPcs <- ncomp
if(ppcaOutput$numIter == opts$maxiters)
{
  print('Maximum number of iterations reached')
}
R2cum      <- rep(NA, nPcs)
TSS        <- sum(myMatsaved^2, na.rm = TRUE)

# for (i in 1:nPcs) {
#   difference <- myMatsaved - (ppcaOutput$scores[,1:i, drop=FALSE] %*% t(ppcaOutput$W[,1:i, drop=FALSE]) )
#   R2cum[i]   <- 1 - (sum(difference^2, na.rm = TRUE) / TSS)
# }

pcaMethodsRes           <- new("pcaRes")
pcaMethodsRes@scores    <- ppcaOutput$scores 
pcaMethodsRes@loadings  <- ppcaOutput$W
pcaMethodsRes@R2cum     <- rep(1, ncomp)
pcaMethodsRes@method    <- "vbpca"

# Return standard ppcaNet output:

output <- list()
output[["W"]]              <- ppcaOutput$W
output[["sigmaSq"]]        <- ppcaOutput$ss
output[["Sigma"]]          <- ppcaOutput$C
output[["numIter"]]  <- ppcaOutput$numIter
output[["pcaMethodsRes"]]  <- pcaMethodsRes

return(output)
}
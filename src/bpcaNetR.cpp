# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' Bayesian PCA updates
//' 
//' Perform parameter updates for PPCA using the Variational Bayes framework
//' from Oba (2003). Not recommended to use standalone, rather it is called from within
//' \code{\link{bpcapM}} and its wrapper \code{\link{pcapM}}.
//' 
//' @param myMat \code{matrix} -- data matrix with observations in rows and 
//' variables in columns. (Note that this is the transpose of \code{X} in
//' \code{\link{pca_full}}.)
//' @param covy \code{matrix} -- the unbiased sample covariance of the data matrix.
//' @param N \code{numeric} -- the number of observations.
//' @param D \code{numeric} -- the number of variables.
//' @param hidden \code{numeric} -- indices of missing values in \code{1:length(myMat)}.
//' @param numberOfNonNAvaluesInEachCol \code{numeric} -- number of observed values in each 
//' column of the data (i.e. variables).
//' @param nomissIndex \code{numeric} -- indices of rows (observations) without any
//' missing values.
//' @param missIndex \code{numeric} -- indices of rows (observations) with missing
//' values.
//' @param nMissing \code{numeric} -- total number of missing values.
//' @param nPcs \code{numeric} -- number of components/latent variables to use.
//' @param threshold \code{numeric} -- threshold for convergence, applied to the precision
//' parameter \code{tau}. Updates for which the change in \code{tau} are below this threshold
//' value stop the algorithm.
//' @param maxIterations \code{numeric} -- the maximum number of iterations to be completed.
//' @return {A \code{list} of 5 elements:
//' \describe{
//' \item{W}{\code{matrix} -- the estimated loadings.}
//' \item{ss}{\code{numeric} -- the estimated model variance.}
//' \item{C}{\code{matrix} -- the estimated covariance matrix.}
//' \item{scores}{\code{matrix} -- the estimated scores.}
//' \item{m}{\code{numeric} -- the estimated mean vector.}
//' }}
//' @seealso \code{\link{bpcapM}}, \code{\link{pcapM}}
//' @examples
//' set.seed(102)
//'   N <- 20
//'   D <- 20
//'   nPcs <- 2
//'   maxIterations <- 1000
//'   X <- matrix(rnorm(50), D, N)
//'   X <- scale(X, center=TRUE, scale=FALSE) # mean 0
//'   covX <- cov(X)
//'   IX <- sample(1:D, 10)
//'   JX <- sample(1:N, 10)
//'   nMissing <- length(IX)+length(JX)
//'   X[JX, IX] <- 0
//'   hidden <- which(X==0)
//'   numberOfNonNAvaluesInEachCol <- colSums(X!=0)
//'   nomissIndex <- which(rowSums(X!=0)==N)
//'   missIndex <- which(rowSums(X!=0)!=N)
//'   threshold <- 1e-4
//'   bpcaNetOutput <- bpcaNet(myMat=X, covy=covX, N=N, D=D, hidden=hidden,
//'     numberOfNonNAvaluesInEachCol=numberOfNonNAvaluesInEachCol,
//'     nomissIndex=nomissIndex, missIndex=missIndex, nMissing=nMissing,
//'     nPcs=nPcs, threshold=threshold, maxIterations=maxIterations)
//' @references Oba, S., Sato, M.A., Takemasa, I., Monden, M.,
//'  Matsubara, K.I. and Ishii, S., 2003.
//'  \href{https://doi.org/10.1093/bioinformatics/btg287}{doi}.
//'  
//'  Stacklies, W., Redestig, H., Scholz, M., Walther, D. and 
//'  Selbig, J., 2007.
//'  \href{https://doi.org/10.1093/bioinformatics/btm069}{doi}.
// [[Rcpp::export()]]
List bpcaNet (arma::mat myMat,
              arma::mat covy,
              int N,
              int D,
              arma::uvec hidden,
              arma::uvec numberOfNonNAvaluesInEachCol,
              arma::uvec nomissIndex,
              arma::uvec missIndex,
              int nMissing,
              int nPcs=2,
              double threshold=1e-4,
              int maxIterations=200) {
  
  //double ss, ss_old, rel_ch, objective, TSS;
  
  double tau, dtau, tauold, trS, ss;
  int nNotMissing = nomissIndex.n_elem;
  int ind;
  
  arma::mat      dy(nNotMissing, D);
  arma::mat      x(nPcs, nNotMissing);
  arma::mat      xTranspose(nNotMissing, nPcs);
  
  arma::mat      proj(N, nPcs);
  arma::mat      yest(N, D);
  arma::mat      myMatNA = myMat;
 // arma::mat      covy(D, D);
  arma::mat      U(D, nPcs);
  arma::mat      Ufull(D, D);
  arma::mat      Vfull(D, D);
  arma::vec      Sfull(D);
  arma::mat      W(D, nPcs);
  arma::vec      S(nPcs);
  arma::mat      V(D, nPcs);
  arma::mat      Rx(nPcs, nPcs);
  arma::mat      Rxinv(nPcs, nPcs);
  arma::mat      Dw(nPcs, nPcs);
  arma::mat      Dwinv(nPcs, nPcs);
  
  arma::mat      identity_nPcs(nPcs, nPcs, arma::fill::eye);
  arma::rowvec   mu(D);
  arma::rowvec   myMatRow(D);
  arma::vec      alpha(nPcs);
  arma::mat      T(D, nPcs);
  yest = myMat;
  arma::field<arma::uvec> nomissidx(N,1);
  arma::field<arma::uvec> missidx(N,1);
  arma::mat covEst(D ,D);
  arma::mat scores(N, nPcs);
  
  scores.fill(arma::datum::nan);
  //  M$scores <- matrix(NA, M$rows, M$comps)

  
  if(nMissing > 0)
  {
    hidden -= 1; // Since R indexing starts from 1, need to convert to C++ indexing
    myMatNA.elem(hidden).fill(arma::datum::inf);
    for(int i =0; i < N; i++)
    {
      nomissidx(i,0) = arma::find_finite(myMatNA.row(i));
      missidx(i,0)   = arma::find_nonfinite(myMatNA.row(i));
    }
  }
  nomissIndex -= 1;
  missIndex   -= 1;

//  covy = arma::cov(yest);
  
  
  arma::svd(Ufull, Sfull, Vfull, covy);
  
  U   = Ufull.cols(0, nPcs - 1 );
  V   = Vfull.rows(0, nPcs - 1 );
  S   = Sfull.subvec(0, nPcs - 1);
  

  mu  = arma::sum(myMat)/numberOfNonNAvaluesInEachCol.t();
  W   = U*arma::diagmat(arma::sqrt(S));

  tau = 1/( trace(covy) - arma::accu(S) );    
  
  double taumax = 1e10;
  double taumin = 1e-10;
  tau  = std::max( std::min( tau, taumax), taumin );
  
  double galpha0 = 1e-10;
  double balpha0 = 1;

  alpha = (2*galpha0 + D)/(tau*arma::diagvec(W.t()*W)+2*galpha0/balpha0);

  double gmu0  = 0.001;
  
  double btau0 = 1;
  double gtau0 = 1e-10;
  
  arma::mat SigW(nPcs, nPcs, arma::fill::eye);
  tauold = 1000;

  for(int i = 0; i < maxIterations; i++)
  {

    Rx    = identity_nPcs + tau*W.t()*W + SigW;
    Rxinv = arma::inv(Rx);
    //std::cout << "A: " << arma::size(myMat.rows(nomissIndex)) << "\n";
    //std::cout << "B: " << arma::size(arma::repmat(mu, nomissIndex.n_elem, 1)) << " " << nomissIndex.n_elem << "\n";
    dy    = myMat.rows(nomissIndex)- arma::repmat(mu, nomissIndex.n_elem, 1);
    
    x     = tau*Rxinv*W.t()*dy.t();
    T     = dy.t()*x.t();
    trS   = arma::accu(dy%dy);

    //Assign the scores for complete rows
    xTranspose = x.t();

    scores.rows(nomissIndex) = xTranspose;
    //std::cout << scores(0,0) << "\n";

    for(int j = 0; j < missIndex.n_elem; j++)
    {

      ind = missIndex(j);
      arma::uvec noMissInds = nomissidx(ind,0);
      arma::uvec missInds   = missidx(ind,0);
      myMatRow         = myMat.row(ind);

      arma::rowvec dyo = myMatRow(noMissInds).t() - mu(noMissInds).t();
      arma::mat Wm     = W.rows(missInds);
      arma::mat Wo     = W.rows(noMissInds);
      Rxinv            = arma::inv(Rx - tau*Wm.t()*Wm);
      arma::vec ex     = tau*Wo.t()*dyo.t();
      arma::vec xx     = Rxinv*ex;
      arma::vec dym    = Wm*xx;

      myMatRow(noMissInds) = dyo;
      myMatRow(missInds)   = dym.t();
      yest.row(ind) = myMatRow + mu;
      T += myMatRow.t()*xx.t();
      T.rows(missInds) += Wm*Rxinv;
      trS += arma::dot(myMatRow, myMatRow) + (missInds.n_elem/tau) + arma::trace(Wm*Rxinv*Wm.t());
      // Assign the scores for rows containing missing values
      scores.row(ind) = xx.t();

    }

    T = T/N;
    trS = trS/N;

    Rxinv = arma::inv(Rx);
    Dw    = Rxinv + tau*T.t()*W*Rxinv + arma::diagmat(alpha/N);
    Dwinv = arma::inv(Dw);
    W     = T*Dwinv;
    tau   = (D + 2*gtau0/N)/(trS - arma::trace(T.t()*W) + (arma::dot(mu, mu)*gmu0+2*gtau0/btau0)/N );
    SigW  = Dwinv*((double)D/(double)N);
    alpha = (2*galpha0+D)/(tau*arma::diagvec(W.t()*W) +arma::diagvec(SigW) + 2*galpha0/balpha0);

    if((i+1) % 10 == 0)
    {
      dtau = std::abs(std::log10(tau) - std::log10(tauold));
      if ( dtau < threshold )
      {
        break;
      }
      tauold = tau;
    }
  }

  ss            = 1/tau;
  covEst        = W*W.t() + (ss*arma::eye<arma::mat>(D,D));

  //std::cout << scores(0,0) << " " << scores(0,1) << " " << scores(0,2) << "\n";
  //std::cout << scores(1,0) << " " << scores(1,1) << " " << scores(1,2) << "\n";
  
  List ret ;
  ret["W"]      = W;
  ret["ss"]     = 1/tau;
  ret["C"]      = covEst;
  ret["scores"] = scores;
  ret["mu"] = mu;
  
  return(ret) ;
}

# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export()]]
List bpcaNet (arma::mat myMat, int N, int D, arma::uvec hidden, arma::uvec numberOfNonNAvaluesInEachCol, arma::uvec nomissIndex, arma::uvec missIndex, int nMissing, int nPcs=2, double threshold=1e-4, int maxIterations=200) {
  
  //double ss, ss_old, rel_ch, objective, TSS;
  
  double tau, dtau, tauold, trS, ss;
  int nNotMissing = N - nMissing;
  int ind;
  
  arma::mat      dy(nNotMissing, D);
  arma::mat      x(nPcs, nNotMissing);
  arma::mat      proj(N, nPcs);
  arma::mat      yest(N, D);
  arma::mat      myMatNA = myMat;
  arma::mat      covy(D, D);
  arma::mat      U(D, nPcs);
  arma::mat      Ufull(D, D);
  arma::mat      Vfull(D, D);
  arma::vec      Sfull(D);
  arma::mat      W(D, nPcs);
  //std::cout << D << " " << nPcs << "\n";
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

  covy = arma::cov(yest);
  
  
  arma::svd(Ufull, Sfull, Vfull, covy);
  
  U   = Ufull.cols(0, nPcs - 1 );
  V   = Vfull.rows(0, nPcs - 1 );
  S   = Sfull.subvec(0, nPcs - 1);
  

  mu  = arma::sum(myMat)/numberOfNonNAvaluesInEachCol.t();
  W   = U*arma::diagmat(arma::sqrt(S));
  //std::cout << W.n_rows <<  " " << W.n_cols <<"\n";
  
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
    //std::cout << i << "\n";
    
    Rx    = identity_nPcs + tau*W.t()*W + SigW;
    Rxinv = arma::inv(Rx);
    dy    = myMat.rows(nomissIndex)- arma::repmat(mu, nNotMissing, 1);
    x     = tau*Rxinv*W.t()*dy.t();
    T     = dy.t()*x.t();
    trS   = arma::accu(dy%dy);

    for(int j = 0; j < nMissing; j++)
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
    //std::cout << "tau: "  << tau << "\n";
    
    
    if((i+1) % 10 == 0)
    {
      dtau = std::abs(std::log10(tau) - std::log10(tauold));
      if (dtau < 0.0001)
      {
        break; 
      }
      tauold = tau;
    }
  }
  
  ss            = 1/tau;
  covEst        = W*W.t() + (ss*arma::eye<arma::mat>(D,D));
  //std::cout << W.n_rows <<  " " << W.n_cols <<"\n";
  
  List ret ;
  ret["W"]      = W;
  ret["ss"]     = 1/tau;
  ret["C"]      = covEst;
  
  return(ret) ;
}

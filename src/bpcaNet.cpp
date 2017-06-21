# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


arma::mat mrdivide(arma::mat A, arma::mat B)
{
    return (solve(B.t(), A.t(), arma::solve_opts::fast )).t();
}

// [[Rcpp::export()]]
List bpcaNet (arma::mat myMat, int N, int D, arma::uvec hidden, arma::uvec numberOfNonNAvaluesInEachCol, arma::uvec nomissIndex, arma::uvec missIndex, int nMissing, int nPcs=2, double threshold=1e-4, int maxIterations=200) {
  
  //double ss, ss_old, rel_ch, objective, TSS;
  
  double tau, dtau, tauold, trS;
  int nNotMissing = N - nMissing;
  int ind;
  
  arma::mat      dy(nNotMissing, D);
  arma::mat      x(nPcs, nNotMissing);
  arma::mat      proj(N, nPcs);
  arma::mat      yest(N, D);
  arma::mat      myMatNA = myMat;
  //arma::sp_mat   covy(D, D);
  arma::mat      covy(D, D);
  arma::mat      U(D, nPcs);
  arma::mat      Ufull(D, D);
  arma::mat      Vfull(D, D);
  arma::vec      Sfull(D);
  arma::mat      W(D, nPcs);
  arma::vec      S(nPcs);
  arma::mat      V(D, nPcs);
  arma::mat      Rx(nPcs, nPcs);
  arma::mat      Rxinv(nPcs, nPcs);
  arma::mat      identity_nPcs(nPcs, nPcs, arma::fill::eye);
  arma::rowvec   mu(D);
  arma::rowvec   myMatRow(D);
  arma::vec      alpha(nPcs);
  arma::mat      T(D, nPcs);
  yest = myMat;
  arma::field<arma::uvec> nomissidx(N,1);
  arma::field<arma::uvec> missidx(N,1);
  
  //std::cout << yest << "\n";
  if(nMissing > 0)
  {
    hidden -= 1; // Since R indexing starts from 1, need to convert to C++ indexing
    myMatNA.elem(hidden).fill(arma::datum::inf);
    for(int i =0; i < N; i++)
    {
      nomissidx(i,0) = arma::find_finite(myMatNA.row(i));
      missidx(i,0)   = arma::find_nonfinite(myMatNA.row(i));
      //std::cout << i << "  here:   " << nomissidx(i,0).t() << "\n";
    }
  }
  nomissIndex -= 1;
  missIndex   -= 1;
  //std::cout << hidden << "\n";
  
  //std::cout << yest.elem(hidden) << "\n";
  //yest.elem(hidden).fill(0);  // Already taken care of in wrapper R function
  //yest.elem(hidden) = arma::zeros<arma::vec>(nMissing);
  //std::cout << yest << "\n";
  
  covy = arma::cov(yest);
  
  
  //arma::svds(U, S, V, covy, nPcs);

  arma::svd(Ufull, Sfull, Vfull, covy);
  
  U   = Ufull.cols(0, nPcs - 1 );
  V   = Vfull.rows(0, nPcs - 1 );
  S   = Sfull.subvec(0, nPcs - 1);
  
  //std::cout << U.n_rows << "\n" << U.n_cols << "\n" << S.n_elem << "\n";

  std::cout << S "\n";
  std::cout << U "\n";
  std::cout << V "\n";
  
  mu  = arma::sum(myMat)/numberOfNonNAvaluesInEachCol.t();
  W   = U*arma::diagmat(arma::sqrt(S));
  tau = 1/( trace(covy) - arma::accu(S) );    
  
  double taumax = 1e10;
  double taumin = 1e-10;
  tau  = std::max( std::min( tau, taumax), taumin );
  
  double galpha0 = 1e-10;
  double balpha0 = 1;

  alpha = (2*galpha0 + D)/(tau*arma::diagvec(W.t()*W)+2*galpha0/balpha0);
  //alpha = (2*galpha0 + D)/(tau*arma::diagvec(W.t()*W)+2*galpha0/balpha0);
  
  //double gmu0  = 0.001;
  
  //double btau0 = 1;
  //double gtau0 = 1e-10;
  arma::mat SigW(nPcs, nPcs, arma::fill::eye);
  tauold = 1000;
  
  for(int i = 0; i < maxIterations; i++)
  {
    std::cout << i << "\n";
    
    Rx    = identity_nPcs + tau*W.t()*W + SigW;
    Rxinv = arma::inv(Rx);
    dy    = myMat.rows(nomissIndex)- arma::repmat(mu, nNotMissing, 1);
    x     = tau*Rxinv*W.t()*dy.t();
    T     = dy.t()*x.t();
    trS   = arma::accu(dy%dy);
    //std::cout << trS << "\n";
    for(int j = 0; j < nMissing; j++)
    {
      ind = missIndex(j);
      arma::uvec noMissInds = nomissidx(ind,0);
      arma::uvec missInds   = missidx(ind,0);
      //std::cout << noMissInds.t() << "\n";
      myMatRow         = myMat.row(ind);
      
      //arma::rowvec   dyo(nPcs);
      arma::rowvec dyo = myMatRow(noMissInds).t() - mu(noMissInds).t();
      arma::mat Wm     = W.rows(missInds);
      arma::mat Wo     = W.rows(noMissInds);
      Rxinv            = arma::inv(Rx - tau*Wm.t()*Wm);
      arma::vec ex     = tau*Wo.t()*dyo.t();
      arma::vec xx     = Rxinv*ex;
      arma::vec dym    = Wm*xx;
      myMatRow(noMissInds) = dyo;
      myMatRow(missInds)   = dym.t();
      
//      std::cout << ex.t() << "\n";
//      std::cout << xx.t() << "\n";
//      std::cout << dym.t() << "\n";
      //std::cout << myMatRow << "\n"<< "\n";
      
      
      
      //arma::rowvec dyoFull = myMat.cols(noMissInds).row(ind);
      //dyo = myMat(ind, nomissidx(j,0));//myMat.row(ind).cols(nomissidx(j,0));
        //myMat(ind, nomissidx(j)) - mu(nomissidx(j));
      //std::cout << myMatRow << "\n";
      
    }
      
    
    if(i == 0){
      //std::cout << nomissIndex << "\n";
      //std::cout << Rx << "\n";
      //std::cout << Rxinv << "\n";
      //std::cout << trS << "\n";
    }
    
    if((i+1) % 10 == 0)
    {
      dtau = std::abs(std::log10(tau) - std::log10(tauold));
      //std::cout << i << "   " << dtau << "\n";
      if (dtau < 0.0001)
      {
        break; 
      }
      tauold = tau;
    }
  }
  
  //std::cout << tau   << "\n";
  //std::cout << alpha << "\n";
  
  
  List ret ;
  ret["W"] = W;
  
  return(ret) ;
}

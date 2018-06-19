# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


arma::mat mrdivide(arma::mat A, arma::mat B)
{
    return (solve(B.t(), A.t(), arma::solve_opts::fast )).t();
}

// [[Rcpp::export()]]
List ppcaNet (arma::mat myMat, int N, int D, arma::mat W, arma::uvec hidden, int nMissing, int nPcs=2, double threshold=1e-5, int maxIterations=1000) {
  
  double ss, ss_old, rel_ch, objective;
  arma::mat    WtW(nPcs, nPcs);
  arma::mat      M(nPcs, nPcs);
  arma::mat      U(nPcs, nPcs);
  arma::mat      V(nPcs, nPcs);
  arma::mat     Sx(nPcs, nPcs);
  arma::mat SumXtX(nPcs, nPcs);
  arma::mat      X(N, nPcs);
  arma::mat   proj(N, nPcs);
  arma::mat  recon(N, D);
  arma::mat covEst(D ,D);

  
  if(nMissing > 0)
  {
    hidden -= 1; // Since R indexing starts from 1, need to convert to C++ indexing
  }
//  std::cout << hidden << "\n";
//  std::cout << myMat(hidden) << "\n";
  
  WtW           = W.t() * W;
  X             = myMat * mrdivide(W, WtW); 
  recon         = X* W.t();
  recon.elem(hidden).fill(0);
  ss = arma::accu(arma::square(recon - myMat)) / (N * D - nMissing); // initial guess for sigma^2
  
  int count = 1, numberOfIterations;
  double old = arma::datum::inf;
  while (count > 0) {
    M  = WtW + ss*arma::eye(nPcs, nPcs);
    U  = arma::chol(M);
    V  = arma::inv(U);
    
    Sx = V*V.t()*ss;  //This is Sigma = inv(I + CC^T/sigma^2) 
    
    //E-step, (co)variances
    ss_old = ss;
    if(nMissing > 0) {
      proj = X * W.t();//See Equation 8 of ppca_mv.pdf
      myMat(hidden) = proj(hidden);//This is y_bar in ppca_mv.pdf
    }
    // E step: expected values
    X = myMat * W * Sx / ss;
    
    //M-step
    SumXtX = X.t() * X;
    W      = mrdivide(myMat.t()*X, (SumXtX + N*Sx)); // Equation 13 of ppca_mv.pdf
    WtW    = W.t()*W;                                
    ss     = ( arma::accu( arma::square(W*X.t() - myMat.t()) ) + N * arma::accu(WtW * Sx) + nMissing * ss_old ) / (N * D);//Equation 14 of ppca_mv.pdf

    
    
    // Note: in the below, we make use of log(det(Sx)) = 2*log(det(V)) - nPcs*log(ss_old)
    objective = N*(D*log(ss)+arma::trace(Sx)-(nPcs*log(ss_old) + 2*log(arma::det(V))) ) +arma::trace(SumXtX) - nMissing*log(ss_old);
    //objective <- N * (D * log(ss) + sum(diag(Sx)) - log(det(Sx)) ) + sum(diag(SumXtX)) - nMissing * log(ss_old)
      
    rel_ch = std::abs( 1 - objective / old );
    old    = objective;
    
    numberOfIterations = count;  
    count += 1;
    if( rel_ch < threshold & count > 5 ) {
      count = 0;
    }
    else if (count > maxIterations) {
      count = 0;
    }
    
  }
  
  covEst        = W*W.t() + (ss*arma::eye<arma::mat>(D,D));
  
  // //log-likelihood
  // double logLikeObs=0;
  // double logLikeProj=0;
  // for (int i=0; i<N; i++){
  //   // calculate the log-likelihood for the observed data
  //   myMat.row(i);
  //   
  //   // calculate the log-likelihood for the projected
  //   
  // }
  
  List ret ;
  ret["W"]      = W;
  ret["ss"]     = ss;
  ret["C"]      = covEst;
  ret["myMat"]  = myMat;
  
  // For us, we don't strictly need the below, because we are only interested in estimating the covariance matrix
  // Leave it in for now (but commented), in case we later decide to give this functionality to the user
  /*
  double TSS;
  arma::mat   differenceMat(N,D);
  arma::mat   eigenVectors(nPcs, nPcs);
  arma::vec   eigenValues(nPcs);
  arma::vec   R2cum(nPcs);
  W = arma::orth(W);

  arma::eig_sym(eigenValues, eigenVectors, arma::cov(myMat * W));
  
  eigenVectors = arma::fliplr(eigenVectors);
  std::reverse(eigenValues.begin(), eigenValues.end());
  
  std::cout << eigenValues  << "\n";
  std::cout << eigenVectors << "\n";
  TSS = 0;
  
  W = W * eigenVectors;
  X = myMat *W;
  
   
  

  TSS   = arma::accu(arma::square(myMat));
  
  for (int i = 0; i < nPcs; i++) {
    differenceMat = myMat - (X.cols(0,i) * arma::trans(W.cols(0,i)));
    R2cum[i] = 1 - (arma::accu(arma::square(differenceMat)) / TSS);
  }
  
  // returns
  ret["scores"]             = X;
  ret["loadings"]           = W;
  ret["myMat"]              = myMat;
  ret["recon"]              = recon;
  ret["proj"]               = proj;
  ret["objective"]          = objective;
  ret["numberOfIterations"] = numberOfIterations;
  ret["R2cum"]              = R2cum;*/
  return(ret) ;
}

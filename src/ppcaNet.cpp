# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


arma::mat mrdivide(arma::mat A, arma::mat B)
{
    return (solve(B.t(), A.t(), arma::solve_opts::fast )).t();
}

// [[Rcpp::export()]]
List ppcaNet (arma::mat myMat, int N, int D, arma::mat C, arma::uvec hidden, int nMissing, int nPcs=2, double threshold=1e-5, int maxIterations=1000) {
  
  double ss, ss_old, rel_ch, objective, TSS;
  arma::mat    CtC(nPcs, nPcs);
  arma::mat      M(nPcs, nPcs);
  arma::mat      U(nPcs, nPcs);
  arma::mat      V(nPcs, nPcs);
  arma::mat     Sx(nPcs, nPcs);
  arma::mat SumXtX(nPcs, nPcs);
  arma::mat      X(N, nPcs);
  arma::mat   proj(N, nPcs);
  arma::mat  recon(N, D);

  arma::mat differenceMat(N,D);
  arma::mat  eigenVectors(nPcs, nPcs);
  arma::vec   eigenValues(nPcs);
  arma::vec   R2cum(nPcs);
  
  if(nMissing > 0)
  {
    hidden -= 1; // Since R indexing starts from 1, need to convert to C++ indexing
  }
//  std::cout << hidden << "\n";
//  std::cout << myMat(hidden) << "\n";
  
  CtC           = C.t() * C;
  X             = myMat * mrdivide(C, CtC); 
  recon         = X* C.t();
  recon.elem(hidden).fill(0);
  ss = arma::accu(arma::square(recon - myMat)) / (N * D - nMissing);
  
  int count = 1, numberOfIterations;
  double old = arma::datum::inf;
  while (count > 0) {
    M  = CtC + ss*arma::eye(nPcs, nPcs);
    U  = arma::chol(M);
    V  = arma::inv(U);
    
    Sx = V*V.t()*ss;
    
    //E-step, (co)variances
    ss_old = ss;
    if(nMissing > 0) {
      proj = X * C.t();
      myMat(hidden) = proj(hidden);
    }
    // E step: expected values
    X = myMat * C * Sx / ss;
    
    //M-step
    SumXtX = X.t() * X;
    C      = mrdivide(myMat.t()*X, (SumXtX + N*Sx));
    CtC    = C.t()*C;
    ss     = ( arma::accu( arma::square(C*X.t() - myMat.t()) ) + N * arma::accu(CtC * Sx) + nMissing * ss_old ) / (N * D);

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
  List ret ;
  ret["C"]      = C;
  ret["ss"]     = ss;
  C = arma::orth(C);

  std::cout << arma::eig_sym(arma::cov(myMat * C)) << "\n";
  
  arma::eig_sym(eigenValues, eigenVectors, arma::cov(myMat * C));
  
  eigenVectors = arma::fliplr(eigenVectors);
  std::reverse(eigenValues.begin(), eigenValues.end());
  
  std::cout << eigenValues  << "\n";
  std::cout << eigenVectors << "\n";
  TSS = 0;
  
  C = C * eigenVectors;
  X = myMat *C;
  
   
  

  TSS   = arma::accu(arma::square(myMat));
  
  for (int i = 0; i < nPcs; i++) {
    //std::cout << X.cols(0,i) << "\n";
    //std::cout << arma::trans(C.cols(0,i)) << "\n";
    
    differenceMat = myMat - (X.cols(0,i) * arma::trans(C.cols(0,i)));
    R2cum[i] = 1 - (arma::accu(arma::square(differenceMat)) / TSS);
  }
  
  // returns
  ret["scores"]             = X;
  ret["loadings"]           = C;
  ret["myMat"]              = myMat;
  ret["recon"]              = recon;
  ret["proj"]               = proj;
  ret["objective"]          = objective;
  ret["numberOfIterations"] = numberOfIterations;
  ret["R2cum"]              = R2cum;
  return(ret) ;
}

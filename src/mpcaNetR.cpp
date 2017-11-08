# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


// [[Rcpp::export()]]
List mpcaNet (const arma::mat Y, arma::mat W, arma::uvec hidden, int nMissing, double v, const double traceS, const int MaxIter = 1000, const double TolFun = 1e-6, const double TolX = 1e-6) {
  
  
  // inputs
  const int n = Y.n_cols;
  const int p = Y.n_rows;
  const int k = W.n_cols;
  
  arma::vec y(p); 
  arma::vec nObs(p);
  
  arma::mat covEst(p,p);
  arma::mat  diagv(k,k);
  arma::mat Sigmax(k,k);
  arma::mat      X(k,n);
  arma::mat     WX(p,n);
  arma::mat      M(k,k);
  arma::mat   sumC(k,k);
  arma::mat   Wnew(p,k);
  
  arma::cube C(k,k,n);
  
  arma::uword currentIndex;
  
  arma::vec x(k);
  arma::vec mu(p);
  arma::vec ww(k);
  
  arma::rowvec Yrow(p);
  
  mu.zeros();
  X.zeros();
  C.zeros();
  
  double vnew, nObsTotal, nloglk, nloglk_new, dw;
  // Y  is p x n
  // Y' is n x p
  // W  is p x k
  
  nloglk = arma::datum::inf;
  if(nMissing > 0)
  {
    hidden -= 1; // Since R indexing starts from 1, need to convert to C++ indexing
  }
  
  
  for (int j = 0; j < p; j++){
    arma::uvec indices = arma::find(Y.row(j)!=0);
    nObs(j) = indices.n_elem;   
  }
  nObsTotal = arma::accu(nObs);
  //std::cout << nObs << "\n";
  
  /*
  std::cout << n << " " << k << "\n";
  std::cout << Y << "\n";
  std::cout << hidden << "\n";
  std::cout << Y(hidden) << "\n";
  
  std::cout << W << "\n";*/
  int    iter = 0;;
  
  
  while(iter < MaxIter)
  {
    
    // PC updates
    diagv = v*arma::eye(k,k);
    for (int i = 0; i < n; i++) {
      y = Y.col(i);
      arma::uvec inds = arma::find(y!= 0);
      arma::mat w = W.rows(inds);
      Sigmax = arma::inv(diagv + (w.t()*w));
      C.slice(i) = Sigmax; // dropping factor of v for brevity
      X.col(i) = Sigmax*w.t()*(y(inds) - mu(inds)); // dropping factor of 1/v for brevity
      //std::cout << inds << "\n";
    }
    /*
    arma::colvec c = arma::sum(Y-W*X,1);
    std::cout << Y << "\n";
    std::cout << W << "\n";
    std::cout << X << "\n";
    std::cout << c << "\n";*/
    
    
    // Update mu
    WX = W*X;
    //std::cout << WX << "\n";
    WX.elem(arma::find(Y == 0)).zeros();
    mu = arma::sum(Y-WX,1);
    //std::cout << mu << "\n";
    for (int j = 0; j < p; j++){
      mu(j) = mu(j)/nObs(j);   
    }
    
    
    // Update W
    for (int j = 0; j < p; j++){
      arma::uvec indices = arma::find(Y.row(j)!=0);
      arma::mat Xcols = X.cols(indices);
      sumC.zeros();
      for(int k = 0; k < indices.n_elem; k++){
        currentIndex = indices(k);
        sumC += C.slice(currentIndex);
      }
      M  = Xcols*Xcols.t() + v*sumC;
      Yrow = Y.row(j);
      ww   = Xcols*(Yrow(indices) - mu(j));
      Wnew.row(j) = (arma::solve(M,ww, arma::solve_opts::fast)).t(); // solves Wnew = inv(M)*ww
      //std::cout << Wnew << "\n";
    }
    
    //Update v
    vnew = 0;
    //std::cout << v << "\n";
    for (int i = 0; i < n; i++) {
      y = Y.col(i);
      arma::uvec inds = arma::find(y!= 0);
      arma::mat  Wred = Wnew.rows(inds);
      arma::mat tmpMat = Wred*C.slice(i)*Wred.t();
      vnew = vnew + arma::accu(arma::square(y(inds) - Wred*X.col(i) - mu(inds)) + v*tmpMat.diag() );
    }
    vnew/=nObsTotal;
    
    nloglk_new = 0;
    for(int i = 0; i < n; i++)
    {
      y = Y.col(i);
      arma::uvec inds = arma::find(y!= 0);
      arma::vec centred = y(inds) - mu(inds);
      arma::mat  Wred = Wnew.rows(inds);
      arma::mat Cy = Wred*Wred.t() + vnew*arma::eye<arma::mat>(inds.n_elem,inds.n_elem);
      nloglk_new += (inds.n_elem*log(2*arma::datum::pi) + log(arma::det(Cy)) + arma::trace(arma::solve(Cy,centred*centred.t(), arma::solve_opts::fast)  ))/2;
    }
    //std::cout << nloglk_new << "\n";
    
    dw = max(max(abs(W-Wnew) / (sqrt(arma::datum::eps)+max(max(abs(Wnew))))));
    //std::cout << dw << "\n";
    
    W = Wnew;
    v = vnew;
    
    
    if(dw < TolX)
    {
      break;
    }
    else if( (nloglk - nloglk_new) < TolFun )
    {
      break;
    }
    
    
    nloglk = nloglk_new;
    iter++;
  }

  
  covEst        = Wnew*Wnew.t() + (v*arma::eye<arma::mat>(p,p));
  
  // returns
  List ret ;
  ret["W"]       = Wnew ;
  ret["ss"]       = vnew ;
  ret["C"]      = covEst;
  ret["numIter"] = iter ;
  
  return(ret) ;
}

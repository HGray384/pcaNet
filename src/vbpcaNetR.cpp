# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


// [[Rcpp::export()]]
List vbpcaNet (const arma::mat Y, arma::mat W, arma::uvec hidden, int nMissing, double v, const double traceS, const int MaxIter = 1000, const double TolFun = 1e-6, const double TolX = 1e-6) {
  
  
  // inputs
  const int n = Y.n_cols;
  const int p = Y.n_rows;
  const int k = W.n_cols;
  
  arma::vec y(p); 
  arma::vec nObs(p);
  
  arma::mat    covEst(p,p);
  arma::mat     diagv(k,k);
  arma::mat diagvwinv(k,k);
  arma::mat    Sigmax(k,k);
  arma::mat    Sigmaw(k,k);
  arma::mat         X(k,n);
  arma::mat        WX(p,n);
  arma::mat         M(k,k);
  arma::mat      sumC(k,k);
  arma::mat      sumD(k,k);
  arma::mat      Wnew(p,k);
  
  arma::cube C(k,k,n);
  arma::cube D(k,k,p);
  
  arma::uword currentIndex;
  
  arma::vec x(k);
  arma::vec mbar(p);
  arma::vec mtilde(p);
  arma::vec ww(k);
  arma::vec vw(k);
  arma::vec vwnew(k);
  
  arma::rowvec Yrow(p);
  
  double vnew, vm, vmnew, nObsTotal, nloglk, nloglk_new, dw;
  
  
  // parameter initialisation
  mbar.zeros();
  mtilde.ones(); // initialised as in Ilin and Raiko code
  X.zeros();
  C.zeros();
  vw=1000*arma::ones(k,1); // initialised as in Ilin and Raiko code
  vm=1000; // initialised as in Ilin and Raiko code
  for (int j = 0; j < p; j++) {
    D.slice(j) = arma::eye(k,k); // each slice is Sigmaw for j-th observation
  }
  
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
  
  
  //std::cout << n << " " << k << "\n";
  //std::cout << Y << "\n";
  //std::cout << hidden << "\n";
  //std::cout << Y(hidden) << "\n";
  
  //std::cout << W << "\n";
  int    iter = 0;
  
  
  while(iter < MaxIter)
  {
  
  // update PCs
  diagv = v*arma::eye(k,k);
  for (int i = 0; i < n; i++) {
  y = Y.col(i);
  arma::uvec inds = arma::find(y!= 0);
  arma::mat w = W.rows(inds);
  sumD.zeros(); // from other
  for(int k = 0; k < inds.n_elem; k++){
    currentIndex = inds(k);
    sumD += D.slice(currentIndex);
  }
  Sigmax = arma::inv(diagv + (w.t()*w+v*sumD));
  C.slice(i) = Sigmax; // dropping a factor of v from Ilin and Raiko
  X.col(i) = Sigmax*w.t()*(y(inds) - mbar(inds));
  //std::cout << inds << "\n";
  }
  
  //arma::colvec c = arma::sum(Y-W*X,1);
  //std::cout << Y << "\n";
  //std::cout << W << "\n";
  //std::cout << X << "\n";
  //std::cout << c << "\n";
  
  
  // Update mbar
  WX = W*X;
  //std::cout << WX << "\n";
  WX.elem(arma::find(Y == 0)).zeros();
  mbar = arma::sum(Y-WX,1);
  //std::cout << "mbar:old \n" << mbar << "\n";
  for (int j = 0; j < p; j++){
  mbar(j) = (vm/(nObs(j)*(vm+(v/nObs(j)))))*mbar(j);   
  }
  //std::cout << "mbar:new \n" << mbar << "\n";
  
  // update mtilde
  //std::cout << "mtilde:old \n" << mtilde << "\n";
  for (int j = 0; j < p; j++){
    mtilde(j) = ((v*vm)/(nObs(j)*(vm+(v/nObs(j)))));   
  }
  //std::cout << "mtilde:new \n" << mtilde << "\n";
  
  // update Sigmaw & W
  diagvwinv = arma::inv(arma::diagmat(vw));
  for (int j = 0; j < p; j++) {
    arma::uvec indices = arma::find(Y.row(j)!=0);
    arma::mat Xcols = X.cols(indices);
    sumC.zeros();
    for(int k = 0; k < indices.n_elem; k++){
      currentIndex = indices(k);
      sumC += C.slice(currentIndex);
    }
    //std::cout << " make it here" << "\n";
    M  = Xcols*Xcols.t() + v*sumC;
    Sigmaw = arma::inv(v*diagvwinv + M);
    D.slice(j) = Sigmaw; // dropping a factor of v from Ilin and Raiko
    //std::cout << " sigmaw" << Sigmaw << "\n";
    // now Wnew
    Yrow = Y.row(j);
    ww   = Xcols*(Yrow(indices) - mbar(j));
    Wnew.row(j) = arma::trans(Sigmaw*ww);
    //std::cout << Wnew << "\n";
  }
  //std::cout << "Wnew: " << Wnew << "\n";
  
  //std::cout << " update v" << "\n";
  
  //Update v
  vnew = 0;
  //std::cout << v << "\n";
  for (int i = 0; i < n; i++) {
    y = Y.col(i);
    arma::uvec inds = arma::find(y!= 0);
    arma::mat  Wred = Wnew.rows(inds);
    arma::mat tmpMat = Wred*C.slice(i)*Wred.t();
    for (int j = 0; j < p; j++) {
      arma::mat tmpMat2 = X.col(i).t()*D.slice(j)*X.col(i);
      arma::mat tmpMat3 = C.slice(i)*D.slice(j); // not sure if matrix product is what is wanted
      vnew = vnew + arma::accu(v*tmpMat2.diag() + v*v*arma::trace(tmpMat3));
    }
    vnew = vnew + arma::accu(arma::square(y(inds) - Wred*X.col(i) - mbar(inds)) + mtilde(inds) + v*tmpMat.diag());
  }
  vnew/=nObsTotal;
  //std::cout << "vnew: " << vnew << "\n";
  
  // update vw
  vwnew.zeros();
  for(int l = 0; l < k; l++){
    for(int j = 0; j < p; j++){
      vwnew(l) = vwnew(l) + (Wnew(j, l)*Wnew(j, l))+v*D(l, l, j);
    }
    vwnew(l)/=p;
  }
  //std::cout << "vwnew: " << vwnew << "\n";
  
  // update vm
  vmnew = 0;
  for(int j = 0; j < p; j++){
    vmnew = vmnew + (mbar(j)*mbar(j))+mtilde(j);
  }
  vmnew/=p;
  //std::cout << "vmnew: " << vmnew << "\n";
  
  
  // compute costs for each parameter
  double cost_y;
  double cost_m;
  double cost_W;
  double cost_X;
  
  // cost_y and cost_m are fast
  cost_y = 0.5 * nObsTotal * (1 + log(2*arma::datum::pi*vnew) );
  cost_m = 0.5 * (p * log(vmnew) - arma::accu(arma::log(mtilde)));
  
  // cost W
  double logDetSigmaw = 0;
  for(int j = 0; j < p; j++)
  {
    logDetSigmaw += log(arma::det(D.slice(j)));
  }
  cost_W = 0.5 * (p * arma::accu(arma::log(vwnew)) - logDetSigmaw);
  
  
  // cost X
  double TrSigmax = 0;
  double logDetSigmax = 0;
  for(int i = 0; i < n; i++)
  {
    TrSigmax += arma::trace(C.slice(i));
    logDetSigmax += log(arma::det(C.slice(i)));
  }
  cost_X = 0.5 * (arma::accu(arma::pow(X, 2)) + TrSigmax - logDetSigmax - n*k);
  
  nloglk_new = cost_y + cost_m + cost_W + cost_X;
  
  
  //std::cout << "new cost: " << nloglk_new << "\n";
  
  dw = max(max(abs(W-Wnew) / (sqrt(arma::datum::eps)+max(max(abs(Wnew))))));
  std::cout << "dw: " << dw << "\n";
  
  W = Wnew;
  v = vnew;
  vw = vwnew;
  vm = vmnew;
  
  
  if(dw < TolX)
  {
  break;
  }
  else if( (nloglk - nloglk_new) < TolFun )
  {
  break;
  }
  
  
  nloglk = nloglk_new;
  iter+=1;
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
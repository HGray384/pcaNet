# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


arma::mat mrdivider(arma::mat A, arma::mat B)
{
    return (solve(B.t(), A.t(), arma::solve_opts::fast )).t();
}

// [[Rcpp::export()]]
List ppcaSensible (const arma::mat Y, arma::mat W, double v, const double traceS, const int MaxIter = 1000, const double TolFun = 1e-6, const double TolX = 1e-6) {
    
    // inputs
    const int n = Y.n_cols;
    const int p = Y.n_rows;
    const int k = W.n_cols;
    
    arma::mat covEst(p ,p);
    
    // Y  is p x n
    // Y' is n x p
    // W  is p x k
    
    // containers
    arma::mat    SW(p,k);
    arma::mat  Wnew(p,k);
    arma::mat     M(k,k);
    arma::mat diagv(k,k);
    arma::mat    CC(p,p);
    arma::mat   Xmu(k,n);
    
    double vnew, dw, dv, delta, nloglk_new, nloglk = arma::datum::inf, eps = arma::datum::eps, sqrteps = sqrt(eps), pi = arma::datum::pi;
    int    iter = 0;;
    
    
    while(iter < MaxIter)
    {
        
        
        diagv = v*arma::eye(k,k);
        
        SW    = Y*(Y.t())*W/(n-1); // Y'W
        M     = ((W.t())*W) + diagv;
        Wnew  = mrdivider(SW, diagv + solve(M, W.t()*SW, arma::solve_opts::fast));
        vnew  = (traceS - trace(mrdivider(SW,M)*Wnew.t()))/p;
        
        
        // Check Convergence:
        dw    = max(max(abs(W-Wnew) / (sqrteps+max(max(abs(Wnew))))));
        dv    = std::abs(v-vnew)/(eps+v);
        dw > dv ? delta = dw : delta = dv;
        
        CC    = Wnew*Wnew.t() + (vnew*arma::eye(p,p));
        
        
        W     = Wnew;
        v     = vnew;
        
        
        nloglk_new = (p*log(2*pi) + log(det(CC)) + trace(solve(CC,Y*Y.t(), arma::solve_opts::fast)/(n-1)) )*n/2;
        
        
        if(delta < TolX)
        {
            break;
        }
        else if( (nloglk - nloglk_new) < TolFun )
        {
            break;
        }
        else if ( std::abs(vnew) < sqrteps )
        {
            break;
        }
        
        nloglk = nloglk_new;

        
        iter++;
        
    }

    //Xmu = solve(M,Wnew.t()*Y, arma::solve_opts::fast);
    
    covEst        = Wnew*Wnew.t() + (v*arma::eye<arma::mat>(p,p));
    
    // returns
    List ret ;
    ret["W"]       = Wnew ;
    //ret["Xmu"]     = Xmu;
    ret["ss"]       = vnew ;
    ret["C"]      = covEst;
    
    ret["numIter"] = iter ;
    /*
     ret["dw"]      = dw ;
    ret["nloglk"]  = nloglk_new;
     */
    return(ret) ;
}

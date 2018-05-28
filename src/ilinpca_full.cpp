# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;

// [[Rcpp::export()]]
List pca_updates (arma::mat X,
               double V,
               arma::mat A,
               arma::cube Av,
               arma::vec Va,
               arma::mat S,
               arma::cube Sv,
               arma::vec Mu,
               arma::vec Muv,
               double Vmu,
               const double hpVa,
               const double hpVb,
               const double hpV,
               const int ndata,
               const arma::vec Nobs_i,
               arma::vec Isv,
               arma::mat M,
               const arma::uvec IX,
               const arma::uvec JX,
               double rms,
               arma::mat errMx,
               const int bias = 1,
               const int niter_broadprior =100,
               const int use_prior = 1,
               const int maxiters = 1000) {
  // housekeeping
  const int n = X.n_cols;
  const int p = X.n_rows;
  const int ncomp = A.n_cols;
  arma::vec dMu;
  arma::vec Mu_old;
  arma::vec th;
  arma::mat Mumat(p, n);
  arma::mat zer(p, n);
  zer.fill(0);
  arma::mat A_j(p, ncomp);
  arma::mat Psi(ncomp, ncomp);
  arma::mat S_i(ncomp, n);
  arma::mat Phi(ncomp, ncomp);
  arma::mat Aold(p, ncomp);
  double cost = 0;
  double cost_old = 0;
  double rms_old = rms;
  int nIter;
  while(!Isv.empty()){
    Isv.reset();
  }
  // std::cout << "I make it past housekeeping" << "\n";
  
  for (int iter = 1; iter < maxiters; iter++) {
    // std::cout << "Iteration: " << iter << "\n";
    nIter = iter;
    // % The prior is not updated at the beginning of learning
    // % to avoid killing sources 
    if (use_prior && iter > niter_broadprior){
      // % Update Va, Vmu
      if (bias){
        Vmu = accu( Mu%Mu );
        if (!Muv.empty()){
          Vmu = Vmu + sum(Muv);
        }
        Vmu = (Vmu + 2*hpVa) / (p + 2*hpVb);
      }
      Va = arma::sum( A%A, 0 ).t();
      if (!Av.empty()){
        for (int i = 1; i< p; i++){
          Va = Va + diagvec(Av.slice(i));
        }
      }
      Va = (Va + 2*hpVa) / (p + 2*hpVb);
      // std::cout << "I update VA" << "\n";
    }
    // update Mu
    if (bias){
      dMu = sum(errMx,1) / Nobs_i ;
      if (!Muv.empty()){
        Muv = V / ( Nobs_i + V/Vmu );
        }
        th = 1 / ( 1 + V/(Nobs_i/Vmu) );  
        Mu_old = Mu;
        Mu = th%( Mu + dMu );
        dMu = Mu - Mu_old;
        Mumat = zer.each_col() + Mu;
        X -= Mumat%M;
    }
    // std::cout << "I update Mu" << "\n";
    
    // % Update S
    if (Isv.empty()){
      for (int j = 0; j<n; j++){
        A_j = repmat(M.col(j),1,ncomp) % A;
        Psi = A_j.t() * A_j + V*arma::eye(ncomp, ncomp);
        if (!Av.empty()){
          // int stopind = notMissingI.size();
          // need to get the stop index (last index of notMissingJ)
          // want to iterate over 
          // i = find(M(:,j))'
          // which is the non-zero elements (rows) for a fixed column
          arma::uvec tmpvec = arma::find(M.col(j));
          int stopind = tmpvec.size();
          for (int l = 0; l<stopind; l++){
            int tmpind = tmpvec(l);
            Psi = Psi + Av.slice(tmpind);
          }
        }
        arma::mat invPsi = inv(Psi);
        S.col(j) = invPsi * A_j.t() * X.col(j);
        Sv.slice(j) = V * invPsi;
      }
        // PrintProgress( opts.verbose, j, n2, 'Updating S:' )
      // }else{
      //   for k = 1:nobscomb{
      //     j = obscombj{k}(1);
      //     %A_j = repmat(full(M(:,j)),1,ncomp) .* A;
      //     A_j = repmat(M(:,j),1,ncomp) .* A;
      //     Psi = A_j' * A_j + diag( repmat(V,1,ncomp) );
      //     if ~isempty(Av){
      //       for i = find(M(:,j))'
      //       Psi = Psi + Av{i};
      //     }
      //   }
        // arma::mat invPsi = inv(Psi);
        // Sv{k} = V * invPsi;
        // tmp = invPsi * A_j';
        // for j = obscombj{k}{
        //   S(:,j) = tmp * X(:,j);
        // }
        // %S(:,obscombj{k}) = tmp * X(:,obscombj{k});
        // PrintProgress( opts.verbose, k, nobscomb, 'Updating S:' )
      // }
    }
    // if opts.verbose == 2, fprintf('\r'), end
    
    // if opts.rotate2pca{
    //   [ dMu, A, Av, S, Sv ] = RotateToPCA( ...
    //     A, Av, S, Sv, Isv, obscombj, opts.bias );
    if (bias){
      Mumat = zer.each_col() + Mu;
      X -= Mumat%M;
      Mu = Mu + dMu;
    }
    // }
    // std::cout << "I update S" << "\n";
    
    //% Update A
    // if opts.verbose == 2
    // fprintf('                                              \r')
    //   end
    for (int i = 0; i<p; i++){
      // %S_i = repmat(full(M(i,:)),ncomp,1) .* S;
      S_i = repmat(M.row(i),ncomp,1) % S;
      Phi = S_i * S_i.t() + diagmat(V/Va);
      //j = find(M(i,:))
      arma::uvec tmpvec2 = arma::find(M.row(i));
      int stopind2 = tmpvec2.size();
      // int stopindj = notMissingI.size();
      for (int k = 0 ; k < stopind2; k++){
        int tmpind2 = tmpvec2(k);
        if (Isv.empty()){
          Phi = Phi + Sv.slice(tmpind2);
        }else{
          Phi = Phi + Sv.slice(Isv(tmpind2));
        }
      }
      arma::mat invPhi = inv(Phi);
      A.row(i) = X.row(i) * S_i.t() * invPhi;
      if (!Av.empty()){
        Av.slice(i) = V * invPhi;
      }
      // 
      // PrintProgress( opts.verbose, i, n1, 'Updating A:' )
    }
    // if opts.verbose == 2, fprintf('\r'), end
    //   
    // [rms,errMx] = compute_rms( X, A, S, M, ndata );
    // prms = compute_rms( Xprobe, A, S, Mprobe, nprobe );
    // std::cout << "I update A" << "\n";
    errMx = (X - A*S)%M;
    rms = pow(accu(errMx%errMx)/ndata, 0.5);
    
    //% Update V
    arma::mat sXv(1, 1);
    if (Isv.empty()){
      // std::cout << "I say that Isv is empty" << "\n";
      for (arma::uword r = 0; r<ndata; r++){
        // arma::mat tmpslice = Sv.slice(JX(r));
        // arma::mat tmpmult = tmpslice * A.row(IX(r)).t();
        sXv += A.row(IX(r)) * Sv.slice(JX(r)) * A.row(IX(r)).t();
        if (!Av.empty()){
          sXv += S.col(JX(r)).t() * Av.slice(IX(r)) * S.col(JX(r)) + accu(Sv.slice(JX(r))%Av.slice(IX(r)));
        }
      }
    }else{
      for (int r = 0;r<ndata; r++){
        sXv += A.row(IX(r)) * Sv.slice(Isv(JX(r))) * A.row(IX(r)).t();
        if (!Av.empty()){
          sXv = sXv + S.col(JX(r)).t() * Av.slice(IX(r)) * S.col(JX(r)) + accu(Sv.slice(Isv(JX(r))) % Av.slice(IX(r)));
        }
      }
    }
    // std::cout << "I update sXv" << "\n";
    if (!Muv.empty()){
      sXv = sXv + accu(Muv(IX));
    }
    // std::cout << "I update Muv after sXv" << "\n";
    sXv = sXv + (pow(rms,2.0))*ndata;
    sXv = as_scalar(sXv);
    // std::cout << "sXv=" << sXv << "\n";
    //%V = rms^2 + V/ndata; 
    V = ( as_scalar(sXv) + 2*hpV ) / (ndata + 2*hpV);
    // std::cout << "I update V" << "\n";
    // 
    // t = toc;
    // lc.rms = [ lc.rms rms ]; lc.prms = [ lc.prms prms ];
    // lc.time = [ lc.time t ];
    
    // if ~isempty(opts.cfstop)
    //   cost = ...
    //     cf_full( X, A, S, Mu, V, Av, Sv, Isv, Muv, Va, Vmu,...
    //     M, sXv, ndata );
    // lc.cost = [ lc.cost cost ];
    // end
    
 /////////////////////////////
 // use_prior = all(~isinf(Va));
 // std::cout << "computing cost" << "\n";
 double cost_x = 0;
 double cost_s = 0;
 cost_x = 0.5 / V * as_scalar(sXv) + 0.5*ndata*log( 2*arma::datum::pi*V );
 // std::cout << "cost_x done" << "\n";
 double cost_mu = 0;
 double cost_a = 0;
 
 if (use_prior){
   if (!Muv.empty()){
     cost_mu = 0.5/Vmu*sum(Mu%Mu+Muv) - 0.5*sum(log(Muv)) + p/2*log(Vmu) - p/2;
   }else if (Vmu != 0){
     cost_mu = 0.5/Vmu*sum(Mu%Mu) + p/2*log(2*arma::datum::pi*Vmu);
   }
   // std::cout << "cost_m done" << "\n";
   if (!Av.empty()){
     cost_a = 0.5*sum(sum(A%A,0).t()/Va) + p/2*sum(log(Va)) - p*ncomp/2;
     // std::cout << "cost_a first product" << "\n";
       for (int i = 0; i<p; i++){
         double val; double sign;
         arma::log_det(val, sign, Av.slice(i));
         double ld = val*sign;
         cost_a += 0.5*sum( diagvec(Av.slice(i)) /Va ) - 0.5*ld;
         // std::cout << "cost_a 2nd product" << "\n";
       }
   }else{
     cost_a = 0.5*sum(sum(A%A,1).t()/Va) + p/2*sum(log(2*arma::datum::pi*Va));
   }
   
   
 }else{ 
   // % No prior for A and Mu
   // %    cost = cost - 0.5*sum(sum(log(2*pi*Av))) - n1*ncomp/2;
   if (!Muv.empty()){
     cost_mu = -0.5*sum(log(2*arma::datum::pi*Muv)) - p/2;
   }
   if (!Av.empty()){
     cost_a = -p*ncomp/2*(1+log(2*arma::datum::pi));
     for (int i = 0; i < p ; i++){
       double val; double sign;
       arma::log_det(val, sign, Av.slice(i));
       double ld = val*sign;
       cost_a = cost_a - 0.5*ld;
   }
 }
 
  }
 // std::cout << "cost_m and cost a done" << "\n";
  cost_s = 0.5*accu(S%S);
  if (Isv.empty()){
    for (int j = 0; j< n; j++){
      double val2; double sign2;
      arma::log_det(val2, sign2, Sv.slice(j));
      double ld2 = val2*sign2;
      cost_s = cost_s + 0.5*arma::trace(Sv.slice(j)) - 0.5*ld2;
    }
  }else{
    for (int j = 0; j < n; j++){
      double val2; double sign2;
      arma::log_det(val2, sign2, Sv.slice(Isv(j)));
      double ld2 = val2*sign2;
      cost_s = cost_s + 0.5*arma::trace(Sv.slice(Isv(j))) - 0.5*ld2;
    }
  }
  cost_s -= ncomp*n/2;
  // std::cout << "cost_s done" << "\n";
  cost = cost_mu + cost_a + cost_x + cost_s;
  // std::cout << "Cost: " << cost << " RMS: "<< rms <<"\n";

 ////////////////////////////
    
      
    //   DisplayProgress( dsph, lc )
    // double angleA = subspace(A,Aold);
    // PrintStep( opts.verbose, lc, angleA )
    //   
    //   convmsg = converg_check( opts, lc, angleA );
    // if ~isempty(convmsg)
    //   if use_prior && iter <= opts.niter_broadprior
    //     % if the prior has never been updated: do nothing
    //       elseif opts.verbose, fprintf( '%s', convmsg ), end
    //       break
    //       end
    /////////////////////////////////////
    // convmsg = '';
    // 
    // if angleA < opts.minangle
    //   convmsg = ...
    //     sprintf( [ 'Convergence achieved (angle between subspaces'...
    //     ' smaller than %.2e)\n' ], opts.minangle );
    
    // elseif opts.earlystop && lc.prms(end) > lc.prms(end-1)
    //   convmsg = sprintf( 'Early stopping.\n' );
    
    // elseif isfield( opts, 'rmsstop' ) && ~isempty(opts.rmsstop) ...
    //   && length(lc.rms)-1 > opts.rmsstop(1)
      
    //   numiter = opts.rmsstop(1);
    // abs_tol = opts.rmsstop(2);
    // rel_tol = [];
    // if length(opts.rmsstop) > 2
    // rel_tol = opts.rmsstop(3);
    // end
    //   rms1 = lc.rms(end-numiter);
    // rms2 = lc.rms(end);
    
    // if abs(rms1-rms2) < abs_tol || ...
    //   ( length(opts.rmsstop) > 2 && abs( (rms1-rms2)/rms2 ) < rel_tol )
    //   
    //   convmsg = ...
    //     sprintf( 'Stop: RMS does not change much for %d iterations.\n',...
    //     numiter );
    // end
      
    //   elseif isfield( opts, 'cfstop' ) && ~isempty(opts.cfstop) && ...
    //     length(lc.cost)-1 > opts.cfstop(1)
    //   
    //   numiter = opts.cfstop(1);
    // abs_tol = opts.cfstop(2);
    // rel_tol = [];
    // if length(opts.cfstop) > 2
    // rel_tol = opts.cfstop(3);
    // end
    //   cost1 = lc.cost(end-numiter);
    // cost2 = lc.cost(end);
    if (fabs(rms_old-rms) < 10e-6 && iter > 2){
      // std::cout << "difference in RMS: " << fabs(rms-rms_old) << "\n";
      // std::cout << "RMS function didn't change: exiting..." << "\n";
      // //   convmsg = ...
      //     sprintf( 'Stop: Cost does not change much for %d iterations.\n',...
      //     numiter );
      // end
      //   elseif nargin >4 && sd_iter == 40
      // convmsg = ...
      //   sprintf(...
      //   [ 'Slowing-down stop. ' ...
      //   'You may continue by changing the gradient type.\n' ] );
      break;
    }
    
    if (fabs(cost_old-cost) < 10e-6 && iter > 2){
    //   std::cout << "difference in cost function: " << fabs(cost-cost_old) << "\n";
    //   std::cout << "cost function didn't change: exiting..." << "\n";
    // //   convmsg = ...
    //     sprintf( 'Stop: Cost does not change much for %d iterations.\n',...
    //     numiter );
    // end
    //   elseif nargin >4 && sd_iter == 40
    // convmsg = ...
    //   sprintf(...
    //   [ 'Slowing-down stop. ' ...
    //   'You may continue by changing the gradient type.\n' ] );
    break;
    }
    /////////////////////////////////////
    
    Aold = A;
    cost_old = cost;
    rms_old = rms;
    sXv.zeros();
        // time = clock;
        // if etime(time,time_autosave) > opts.autosave
        //   time_autosave = time;
        // if opts.verbose == 2, fprintf('Saving ... '), end
        //   save( opts.filename,...
        //     'A', 'S', 'Mu', 'V', 'Av', 'Muv', 'Sv', 'Isv',...
        //     'Va', 'Vmu',...
        //     'lc', 'Ir', 'Ic', 'n1x', 'n2x', 'n1', 'n2' )
        //   if opts.verbose == 2, fprintf('done\n'), end
        //     end
        //     end
        // std::cout << "Mu " <<  << "\n";
}
  // % Finally rotate to the PCA solution
  // if ~opts.rotate2pca
  //   [ dMu, A, Av, S, Sv ] = RotateToPCA( ...
  //     A, Av, S, Sv, Isv, obscombj, opts.bias );
  if (bias){ 
    Mu = Mu + dMu;
  }
  // end
  //   end
  
  // if n1 < n1x
  //   [ A, Av ] = addmrows( A, Av, Ir, n1x, Va );
  // [ Mu, Muv ] = addmrows( Mu, Muv, Ir, n1x, Vmu );
  // end
  //   if n2 < n2x
  //     [ S, Sv, Isv ] = addmcols( S, Sv, Ic, n2x, Isv );
  //   end
  
  // A.S = S;
  // A.Mu = Mu;
  // A.V = V;
  // A.Va = Va;
  // A.Vmu = Vmu;
  // A.Av = Av;
  // A.Sv = Sv;
  // A.Isv = Isv;
  // A.Muv = Muv;
  // A.lc = lc;
  arma::mat C = A*A.t() + V*arma::eye(p, p);
  // returns
  List ret;
  ret["scores"] = S.t();
  // ret["m"] = Mu;
  // ret["vm"] = Vmu;
  ret["ss"] = V;
  ret["W"] = A;
  // ret["vw"] = Va;
  // ret["sigmaw"] = Av;
  // ret["sigmascores"] = Sv;
  // ret["Isv"]= Isv;
  // ret["mtilde"] = Muv;
  ret["C"] = C;
  // ret["rms"] = rms;
  ret["numIter"] = nIter;
  
  return(ret) ;
}
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
                  const int rotate2pca = 1,
                  const int niter_broadprior =100,
                  const int use_prior = 1,
                  const int use_postvar = 1,
                  const int maxiters = 1000,
                  const int verbose = 1) {
  // housekeeping
  const int n = X.n_cols;
  const int p = X.n_rows;
  const int ncomp = A.n_cols;
  arma::vec dMu;
  arma::vec Mu_old;
  arma::vec th;
  arma::vec lcrms(maxiters + 1);
  arma::vec lccost(maxiters + 1);
  arma::mat Mumat(p, n);
  arma::mat zer(p, n);
  zer.fill(0);
  arma::mat A_j(p, ncomp);
  arma::mat Psi(ncomp, ncomp);
  arma::mat S_i(ncomp, n);
  arma::mat Phi(ncomp, ncomp);
  arma::mat Aold(p, ncomp);
  double cost = 0;
  double cost1, cost2, cost_old = 0;
  double rms1, rms2, rms_old = rms;
  int nIter;
  int rmsStopNumIter = 100;
  int costStopNumIter = 100;
  Aold = A;
  lcrms(0)  = rms;
  lccost(0) = arma::datum::nan;
  while(!Isv.empty()){
    Isv.reset();
  }
  //std::cout << "I make it past housekeeping" << "\n";

  if(use_postvar == 0)
  {
    Muv.clear();
    Av.clear();
  }
  
  if (bias == 0)
  {
    Muv.clear();
  }
  
  for (int iter = 1; iter < maxiters; iter++) {
    // std::cout << "Iteration: " << iter << "\n";
    nIter = iter;
    // The prior is not initially updated
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
      // std::cout << "Va: " << Va << "\n";
    }
    // update Mu
    if (bias){
      dMu = sum(errMx,1) / Nobs_i ;
      // Rcout<<"\n dMu at init is : " << std::endl;
      // for(int i=0; i < p; ++i){
      //     Rcout << dMu(i) << " ";
      // }
      // std::cout << "dMu[0]: " << dMu[0] << "\n";
      // std::cout << "dMu[1]: " << dMu[1] << "\n";
      
      if (!Muv.empty()){
        Muv = V / ( Nobs_i + V/Vmu );
        // Rcout<<"\n Muv after Muv if is : " << std::endl;
        // for(int i=0; i < p; ++i){
        //   Rcout << Muv(i) << " ";
        // }
        // std::cout << "Muv[0]: "  << Muv[0] << "\n";
        // std::cout << "Muv[1]: "  << Muv[1] << "\n";
        
      }
      th = 1 / ( 1 + (V/Nobs_i)/Vmu );
      // Rcout<<"\n th is : " << std::endl;
      // for(int i=0; i < p; ++i){
      //   Rcout << th(i) << " ";
      // }
  
      //std::cout << "th[0]: "  << th[0] << "\n";
      //std::cout << "th[1]: "  << th[1] << "\n";
      
      Mu_old = Mu;
      Mu = th%( Mu + dMu );
      // Rcout<<"\n Mu after update is : " << std::endl;
      // for(int i=0; i < p; ++i){
      //   Rcout << Mu(i) << " ";
      // }
      dMu = Mu - Mu_old;
      // Rcout<<"\n new dMu is : " << std::endl;
      // for(int i=0; i < p; ++i){
      //   Rcout << dMu(i) << " ";
      // }
      Mumat = zer.each_col() + dMu;
      // Rcout <<"\n Mumat is : " << std::endl;
      // for(int i=0; i < p; ++i) {
      //   for(int  j=0; j < n; ++j)
      //     Rcout << Mumat(i, j) << " ";
      // }
      
      // std::cout << "accu(X) : "  << accu(X)  << "\n";
      // std::cout << "X[0,0] - Mu[0] : "  << X(0,0) - Mu[0] << "\n";
      X -= Mumat%M;
      //std::cout << "Mumat%M : "  << Mumat%M  << "\n";
      // std::cout << "X[0,0] : "  << X(0,0)  << "\n";
      // std::cout << "X[0,1] : "  << X(0,1)  << "\n";
      // std::cout << "X[1,0] : "  << X(1,0)  << "\n";
      // std::cout << "accu(X) : "  << accu(X)  << "\n";
      
      // Rcout <<"\n X (post mu subtract) is : " << std::endl;
      // for(int i=0; i < p; ++i) {
      //   for(int  j=0; j < n; ++j)
      //     Rcout << X(i, j) << " ";
      // }
    }
    //std::cout << "Mu[0]: "  << Mu[0] << "\n";
    //std::cout << "Mu[1]: "  << Mu[1] << "\n";
    // std::cout << "dMu[0]: " << dMu[0] << "\n";
    // std::cout << "dMu[1]: " << dMu[1] << "\n";
    
    // std::cout << "I update Mu" << "\n";
    
    // % Update S
    if (Isv.empty()){
      for (int j = 0; j<n; j++){
        A_j = repmat(M.col(j),1,ncomp) % A;
        // std::cout << "A_j[0, 0]: "  << A_j(0, 0) << "\n";
        // std::cout << "A_j[0, 1]: "  << A_j(0, 1) << "\n";
        //std::cout << accu(A_j)  << "\n";
        Psi = A_j.t() * A_j + V*arma::eye(ncomp, ncomp);
        //std::cout << " A_j.t() * A_j: "  <<  A_j.t() * A_j << "\n";
        //std::cout << "Psi: "  << Psi << "\n";
        if (!Av.empty()){
          // int stopind = notMissingI.size();
          // need to get the stop index (last index of notMissingJ)
          // want to iterate over 
          // i = find(M(:,j))'
          // which is the non-zero elements (rows) for a fixed column
          arma::uvec tmpvec = arma::find(M.col(j));
          //std::cout << "tmpvec: "  << tmpvec << "\n";
          
          int stopind = tmpvec.size();
          for (int l = 0; l<stopind; l++){
            int tmpind = tmpvec(l);
            Psi = Psi + Av.slice(tmpind);
            //std::cout << "Av.slice(tmpind): "  << Av.slice(tmpind) << "\n";
          }
        }
        //std::cout << "Psi: "  << Psi << "\n";
        arma::mat invPsi = inv(Psi);
        // if( j== 0)
        // {
        //   std::cout << "accu(X.col(j)) : "  << accu(X.col(j))  << "\n";
        //   std::cout << "accu(A_j.t()): "  << accu(A_j.t())  << "\n";
        //   std::cout << "A_j.t() * X.col(j): "  << A_j.t() * X.col(j)  << "\n";
        // }
        //std::cout << "invPsi * A_j.t() * X.col(j): "  << invPsi * A_j.t() * X.col(j) << "\n";
        //std::cout << "invPsi * A_j.t() * X.col(j): "  << invPsi * A_j.t() * X.col(j) << "\n";
        //std::cout << "invPsi * A_j.t() * X.col(j): "  << invPsi * A_j.t() * X.col(j) << "\n";
        
        S.col(j) = invPsi * A_j.t() * X.col(j);
        Sv.slice(j) = V * invPsi;
      }
      
      // std::cout << "Psi: "  << Psi << "\n";
      // 
      // std::cout << "S[0, 0]: "  << S(0, 0) << "\n";
      // std::cout << "S[0, 1]: "  << S(0, 1) << "\n";
      // std::cout << "S[1, 0]: "  << S(1, 0) << "\n";
      // std::cout << accu(S) << "\n";
      // std::cout << "Sv.slice(0): "  << Sv.slice(0) << "\n";
      // std::cout << "Sv.slice(1): "  << Sv.slice(1) << "\n";
      
            
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
    // if (bias){
    //   Mumat = zer.each_col() + Mu;
    //   X -= Mumat%M;
    //   Mu = Mu + dMu;
    // }
    // }
    // Rcout<<"\n S (pre-rotate) is : " << std::endl;
    // for(int i=0; i < ncomp; ++i) {
    //   for(int  j=0; j < n; ++j)
    //     Rcout << S(i, j) << " ";
    // }
    if (rotate2pca){
      //   [ dMu, A, Av, S, Sv ] = RotateToPCA( ...
      //     A, Av, S, Sv, Isv, obscombj, opts.bias );
      if (verbose){
        Rcout << "rotating to PCA basis..." << std::endl;
      }
      if (bias){
        //   mS = mean(S,2);
        // Rcout << "subtracting mean of S" << std::endl;
        arma::vec mS = mean(S, 1);
        // dMu = A*mS;
        dMu = A*mS;
        // S = S - repmat(mS,1,n2);
        S -= repmat(mS, 1, n);
        // Rcout << "mean subtracted" << std::endl;
        // else
      }
      //   dMu = 0;
      // end
      //   
      //   covS = S*S';
      arma::mat covS = S*S.t();
      // if isempty(Isv)
      if (Isv.empty()){
        // Rcout << "adding Sv to covS" << std::endl;
        //   for j = 1:n2
        for (int j = 0; j<n; j++){
          //     covS = covS + Sv{j};
          covS += Sv.slice(j);
          // end
        }
      }
      //   else
      //     nobscomb = length(obscombj);
      //   for j = 1:nobscomb
      //     covS = covS + ( length(obscombj{j})*Sv{j} );
      //   end
      //     end
      //     
      //     covS = covS / n2;
      //   %covS = covS / (n2-n1);
      //   [VS,D] = eig(covS);
      // Rcout << "computing eigenvalue decomposition of covS" << std::endl;
      covS = covS / n;
      arma::mat VS;
      arma::vec eigvals;
      eig_sym(eigvals, VS, covS);
      arma::mat D = diagmat(eigvals);
      //   RA = VS*sqrt(D);
      arma::mat RA = VS*sqrt(D);
      //   A = A*RA;
      A = A*RA;
      //   covA = A'*A;
      arma::mat covA = A.t()*A;
      //   if ~isempty(Av)
      if (!Av.empty()){
        // Rcout << "adding Av to covA" << std::endl;
        //     for i = 1:n1
        for (int i = 0; i<p ; i++){
          //       Av{i} = RA'*Av{i}*RA;
          Av.slice(i) = RA.t()*Av.slice(i)*RA;
          //   covA = covA + Av{i};
          covA += Av.slice(i);
          //   end
        }
        //     end
      }
      //     covA = covA / n1;
      //   [VA,DA] = eig(covA);
      // Rcout << "computing eigen decomp of covA" << std::endl;
      covA = covA / p;
      arma::mat VA;
      arma::vec eigvalsA;
      eig_sym(eigvalsA, VA, covA);
      arma::mat DA = diagmat(eigvalsA);
      // I don't think the below 3 lines need to be translated
      // this is because c++ returns the eigenvalues already sorted
      // by largest to smallest
      //   [DA,I] = sort( -diag(DA) );
      //   DA = -DA;
      //   VA = VA(:,I);
      //   A = A*VA;
      arma::uvec I;
      I = arma::sort_index(-arma::diagvec(DA));
      // Rcout<<"\n I is : " << std::endl;
      // for(int j=0; j < ncomp; ++j) {
      //   Rcout << I(j) << " ";
      // }
      arma::vec DAsortvec;
      DAsortvec = sort(-arma::diagvec(DA));
      // Rcout<<"\n DAsortvec is : " << std::endl;
      // for(int j=0; j < ncomp; ++j) {
      //   Rcout << DAsortvec(j) << " ";
      // }
      DA = arma::diagmat(-DAsortvec);
      // Rcout<<"\n DA is : " << std::endl;
      // for(int i=0; i < ncomp; ++i) {
      //   for(int  j=0; j < ncomp; ++j)
      //     Rcout << DA(i, j) << " ";
      // }
      // Rcout<<"\n VA is : " << std::endl;
      // for(int i=0; i < ncomp; ++i) {
      //   for(int  j=0; j < ncomp; ++j)
      //     Rcout << VA(i, j) << " ";
      // }
      arma::mat VAsorted = VA;
      for (int j = 0; j<ncomp; j++){
        int tmpind = I(j);
        VAsorted.col(j) = VA.col(tmpind);
      }
      // Rcout<<"\n VAsorted is : " << std::endl;
      // for(int i=0; i < ncomp; ++i) {
      //   for(int  j=0; j < ncomp; ++j)
      //     Rcout << VAsorted(i, j) << " ";
      // }
      VA = VAsorted;
      A = A*VA;
      //   
      //   if ~isempty(Av)
      if (!Av.empty()){
        // Rcout << "updating Av" << std::endl;
        //     for i = 1:n1
        for (int i = 0; i<p; i++){
          //       Av{i} = VA'*Av{i}*VA;
          Av.slice(i) = VA.t()*Av.slice(i)*VA;
          //   end
        }
        //     end
      }
      //     R = VA'*diag(1./sqrt(diag(D)))*VS';
      // Rcout << "computing rotation matrix R" << std::endl;
      arma::mat R = VA.t()*diagmat(1/sqrt(diagvec(D)))*VS.t();
      //   
      //   S = R*S;
      S = R*S;
      //   for j = 1:length(Sv)
      // Rcout << "updating Sv" << std::endl;
      for (int j = 0; j<n; j++){
        //     Sv{j} = R*Sv{j}*R';
        Sv.slice(j) = R*Sv.slice(j)*R.t();
        //   end
      }
      if (bias){ 
        // Rcout << "updating Mu" << std::endl;
        Mumat = zer.each_col() + dMu;
        X -= Mumat%M;
        Mu = Mu + dMu;
        // end
      }
      //   end
    }
    
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
    
    
    // std::cout << "A[0, 0]: "  << A(0, 0) << "\n";
    // std::cout << "A[0, 1]: "  << A(0, 1) << "\n";
    // std::cout << "A[1, 0]: "  << A(1, 0) << "\n";
    // std::cout << "A[1, 1]: "  << A(1, 1) << "\n";
    // 
    // std::cout << accu(A) << "\n";
    // std::cout << "Av.slice(0): "  << Av.slice(0) << "\n";
    // std::cout << "Av.slice(1): "  << Av.slice(1) << "\n";
    
    // if opts.verbose == 2, fprintf('\r'), end
    //   
    // [rms,errMx] = compute_rms( X, A, S, M, ndata );
    // prms = compute_rms( Xprobe, A, S, Mprobe, nprobe );
    // std::cout << "I update A" << "\n";
    // Rcout<<"\n X pre-compute_rms is : " << std::endl;
    // for(int i=0; i < p; ++i) {
    //   for(int  j=0; j < n; ++j)
    //     Rcout << X(i, j) << " ";
    // }
    // Rcout<<"\n A pre-compute_rms is : " << std::endl;
    // for(int i=0; i < p; ++i) {
    //   for(int  j=0; j < ncomp; ++j)
    //     Rcout << A(i, j) << " ";
    // }
    // Rcout<<"\n S pre-compute rms is : " << std::endl;
    // for(int i=0; i < ncomp; ++i) {
    //   for(int  j=0; j < n; ++j)
    //     Rcout << S(i, j) << " ";
    // }
    // Rcout<<"\n M pre-compute_rms is : " << std::endl;
    // for(int i=0; i < p; ++i) {
    //   for(int  j=0; j < n; ++j)
    //     Rcout << M(i, j) << " ";
    // }
    errMx = (X - A*S)%M;
    // Rcout<<"\n errMx is : " << std::endl;
    // for(int i=0; i < p; ++i) {
    //   for(int  j=0; j < n; ++j)
    //     Rcout << errMx(i, j) << " ";
    // }
    rms = pow(accu(errMx%errMx)/ndata, 0.5);
    // std::cout << accu(errMx) << "\n";
    // std::cout << "rms: "  << rms << "\n";
    
    //% Update V
    arma::mat sXv(1, 1);
    sXv.zeros();
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
    // Rcout << "sXv after big loop is" << sXv << "\n";
    // std::cout << "I update sXv" << "\n";
    if (!Muv.empty()){
      sXv = sXv + accu(Muv(IX));
    }
    // Rcout << "sXv after Muv 'if' is" << sXv << "\n";
    // std::cout << "I update Muv after sXv" << "\n";
    // Rcout << "rms before power raise is" << rms << "\n";
    sXv = sXv + (pow(rms,2.0))*ndata;
    // Rcout << "sXv after power raise is" << sXv << "\n";
    sXv = as_scalar(sXv);
    // Rcout << "sXv after scalar change is" << sXv << "\n";
    // std::cout << "sXv=" << sXv << "\n";
    //%V = rms^2 + V/ndata; 
    V = ( as_scalar(sXv) + 2*hpV ) / (ndata + 2*hpV);
    
    lcrms(iter) = rms;
    
    // std::cout << "V: "  << V << "\n";
    // std::cout << "sXv: "  << sXv << "\n";

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
        // std::cout << "Here a4" << "\n";
        // std::cout << "sum(A%A,0).t()" << sum(A%A,0).t() << "\n";
        // std::cout << "Va" << Va << "\n";
        // std::cout << "0.5*sum(sum(A%A,0).t()/Va) + p/2*sum(log(2*arma::datum::pi*Va))" << 0.5*sum(sum(A%A,0).t()/Va) + p/2*sum(log(2*arma::datum::pi*Va)) << "\n";
        // 
        cost_a = 0.5*sum(sum(A%A,0).t()/Va) + p/2*sum(log(2*arma::datum::pi*Va));

      }
    }
    else{ 
      //std::cout << "Here b" << "\n";
      
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
    lccost(iter) = cost;
    if (verbose){
      Rcout << "Step: " << iter << ", Cost: " << cost << ", RMS: "<< rms << std::endl;
    }
    ////////////////////////////
    
    // arma::mat U1;
    // arma::mat U2;
    // arma::vec s1;
    // arma::mat V1;
    // 
    // svd_econ(U1, s1, V1, A);
    // svd_econ(U2, s1, V1, Aold);
    // if(iter == 1)
    // {
    //   std::cout << "U1(0,0): " << U1(0,0) << "\n";
    //   std::cout << "U1(0,1): " << U1(0,1) << "\n";
    //   std::cout << "U2(0,0): " << U2(0,0) << "\n";
    //   std::cout << "U2(0,1): " << U2(0,1) << "\n";
    // }
    // if(U1.n_cols < U2.n_cols)
    // {
    //   arma::mat U2new = U1;
    //   arma::mat U1new = U2;
    // }
    // else
    // {
    //   arma::mat U2new = U2;
    //   arma::mat U1new = U1;
    // }
    // U2new = U2new - U1new*(U1new.t()*U2new);
    // theta = asin(min(ones(superiorfloat(A,B)),norm(B)));
    
    //    tmpA     = orth(A);
//    tmpAold  = orth(Aold);
    
    // if()
    // {
    //   
    // }
    // 
    // angleA = subspace(A,Aold);
    
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
    // if (fabs(rms_old-rms) < 10e-6 && iter > 2){
    if (iter > rmsStopNumIter){
      rms1 = lcrms(iter-rmsStopNumIter);
      rms2 = lcrms(iter);
      if(fabs(rms2-rms1) < 1e-6)
      {
        if (verbose){
          Rcout << "Stop: RMS has not changed much for 100 iterations." << std::endl;
        }
        break;
      }
      
      
      cost1 = lccost(iter-costStopNumIter);
      cost2 = lccost(iter);
      if(fabs(cost2-cost1) < 1e-4)
      {
        if (verbose){
          Rcout << "Stop: Cost has not changed much for 100 iterations." << std::endl;
        }
        break;
      }
      
    }
    
    // if (fabs(cost_old-cost) < 10e-6 && iter > 2){
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
      // break;
    // }
    /////////////////////////////////////
    
    Aold = A;
    cost_old = cost;
    rms_old = rms;
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
    // Rcout<<"\n A is : " << A(0, 0) << std::endl;
    // Rcout<<"\n A is : " << std::endl;
    // for(int i=0; i < p; ++i) {
    //   for(int  j=0; j < ncomp; ++j)
    //     Rcout << A(i, j) << " ";
    // }
    // Rcout<<"\n S is : " << std::endl;
    // for(int i=0; i < ncomp; ++i) {
    //   for(int  j=0; j < n; ++j)
    //     Rcout << S(i, j) << " ";
    // }
    // Rcout<<"\n Mu is : " << std::endl;
    // for(int i=0; i < p; ++i) {
    //     Rcout << Mu(i) << " ";
    // }
    // Rcout<<"\n V is : " << std::endl;
    // Rcout << V << std::endl;
    // Rcout<<"\n Muv is : " << std::endl;
    // for(int i=0; i < p; ++i) {
    //   Rcout << Muv(i) << " ";
    // }
  }
  
  if (!rotate2pca){
    //   [ dMu, A, Av, S, Sv ] = RotateToPCA( ...
    //     A, Av, S, Sv, Isv, obscombj, opts.bias );
    if (verbose){
      Rcout << "rotating to PCA basis..." << std::endl;
    }
    if (bias){
      //   mS = mean(S,2);
      // Rcout << "subtracting mean of S" << std::endl;
      arma::vec mS = mean(S, 1);
      // dMu = A*mS;
      dMu = A*mS;
      // S = S - repmat(mS,1,n2);
      S -= repmat(mS, 1, n);
      // Rcout << "mean subtracted" << std::endl;
      // else
    }
    //   dMu = 0;
    // end
    //   
    //   covS = S*S';
    arma::mat covS = S*S.t();
    // if isempty(Isv)
    if (Isv.empty()){
      // Rcout << "adding Sv to covS" << std::endl;
      //   for j = 1:n2
      for (int j = 0; j<n; j++){
        //     covS = covS + Sv{j};
        covS += Sv.slice(j);
        // end
      }
    }
    //   else
    //     nobscomb = length(obscombj);
    //   for j = 1:nobscomb
    //     covS = covS + ( length(obscombj{j})*Sv{j} );
    //   end
    //     end
    //     
    //     covS = covS / n2;
    //   %covS = covS / (n2-n1);
    //   [VS,D] = eig(covS);
    // Rcout << "computing eigenvalue decomposition of covS" << std::endl;
    covS = covS / n;
    arma::mat VS;
    arma::vec eigvals;
    eig_sym(eigvals, VS, covS);
    arma::mat D = diagmat(eigvals);
    //   RA = VS*sqrt(D);
    arma::mat RA = VS*sqrt(D);
    //   A = A*RA;
    A = A*RA;
    //   covA = A'*A;
    arma::mat covA = A.t()*A;
    //   if ~isempty(Av)
    if (!Av.empty()){
      // Rcout << "adding Av to covA" << std::endl;
      //     for i = 1:n1
      for (int i = 0; i<p ; i++){
        //       Av{i} = RA'*Av{i}*RA;
        Av.slice(i) = RA.t()*Av.slice(i)*RA;
        //   covA = covA + Av{i};
        covA += Av.slice(i);
        //   end
      }
      //     end
    }
    //     covA = covA / n1;
    //   [VA,DA] = eig(covA);
    // Rcout << "computing eigen decomp of covA" << std::endl;
    covA = covA / p;
    arma::mat VA;
    arma::vec eigvalsA;
    eig_sym(eigvalsA, VA, covA);
    arma::mat DA = diagmat(eigvalsA);
    // I don't think the below 3 lines need to be translated
    // this is because c++ returns the eigenvalues already sorted
    // by largest to smallest
    //   [DA,I] = sort( -diag(DA) );
    //   DA = -DA;
    //   VA = VA(:,I);
    //   A = A*VA;
    arma::uvec I;
    I = arma::sort_index(-arma::diagvec(DA));
    // Rcout<<"\n I is : " << std::endl;
    // for(int j=0; j < ncomp; ++j) {
    //   Rcout << I(j) << " ";
    // }
    arma::vec DAsortvec;
    DAsortvec = sort(-arma::diagvec(DA));
    // Rcout<<"\n DAsortvec is : " << std::endl;
    // for(int j=0; j < ncomp; ++j) {
    //   Rcout << DAsortvec(j) << " ";
    // }
    DA = arma::diagmat(-DAsortvec);
    // Rcout<<"\n DA is : " << std::endl;
    // for(int i=0; i < ncomp; ++i) {
    //   for(int  j=0; j < ncomp; ++j)
    //     Rcout << DA(i, j) << " ";
    // }
    // Rcout<<"\n VA is : " << std::endl;
    // for(int i=0; i < ncomp; ++i) {
    //   for(int  j=0; j < ncomp; ++j)
    //     Rcout << VA(i, j) << " ";
    // }
    arma::mat VAsorted = VA;
    for (int j = 0; j<ncomp; j++){
      int tmpind = I(j);
      VAsorted.col(j) = VA.col(tmpind);
    }
    // Rcout<<"\n VAsorted is : " << std::endl;
    // for(int i=0; i < ncomp; ++i) {
    //   for(int  j=0; j < ncomp; ++j)
    //     Rcout << VAsorted(i, j) << " ";
    // }
    VA = VAsorted;
    A = A*VA;
    //   
    //   if ~isempty(Av)
    if (!Av.empty()){
      // Rcout << "updating Av" << std::endl;
      //     for i = 1:n1
      for (int i = 0; i<p; i++){
        //       Av{i} = VA'*Av{i}*VA;
        Av.slice(i) = VA.t()*Av.slice(i)*VA;
        //   end
      }
      //     end
    }
    //     R = VA'*diag(1./sqrt(diag(D)))*VS';
    // Rcout << "computing rotation matrix R" << std::endl;
    arma::mat R = VA.t()*diagmat(1/sqrt(diagvec(D)))*VS.t();
    //   
    //   S = R*S;
    S = R*S;
    //   for j = 1:length(Sv)
    // Rcout << "updating Sv" << std::endl;
    for (int j = 0; j<n; j++){
      //     Sv{j} = R*Sv{j}*R';
      Sv.slice(j) = R*Sv.slice(j)*R.t();
      //   end
    }
    if (bias){ 
      // Rcout << "updating Mu" << std::endl;
      Mu = Mu + dMu;
      // end
    }
    //   end
  }

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
  // covariance matrix
  arma::mat C = A*A.t() + V*arma::eye(p, p);
  
  // double detTerm=0;
  // double distTerm=0;
  // double normTerm=0;
  // double val3=0;
  // double sign3=0;
  // double logLikeObs=0;
  // double logLikeproj=0;
  // // calculate log likelihood
  // for (int i=0; i<n; i++){
  //   
  //   if ()
  //   arma::mat truncC;
  //   arma::log_det(val3, sign3, truncC);
  //   detTerm += val3*sign3;
  // }
  
  // returns
  List ret;
  ret["scores"] = S.t();
  ret["m"] = Mu;
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
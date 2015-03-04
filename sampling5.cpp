// zero-inflated hierarchical poisson regression using Gibbs sampling 

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// normal
// [[Rcpp::export]]
arma::vec armaNormal(const int N) {
    arma::vec x = arma::randn(N,1);
    return x;
}

// multivariate normal
// [[Rcpp::export]]
mat mvrnormArma(int n, vec mu, mat sigma) {
   int ncols = sigma.n_cols;
   mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
mat generateTmat(double alpha){
  int T = 11;
  int i = 0; int j = 0;
  vec v(T);
  mat Tmat = eye(T,T);

  for (i=0; i<T; i++){
    v[i] = exp(-i*alpha);}
  
  for (i=0; i<T; i++){
    for (j=0; j<T; j++){
      if (i == j){Tmat(i,j) = v[0];}
      else if (i < j){Tmat(i,j) = v[j-i];}
      else Tmat(i,j) = v[i-j];};}
  
  //List tempres;
  //tempres["Tmmat"] = Tmat;
  //tempres["v"] = v;
  //return tempres;
  return Tmat;}
  
  
// mustar
// [[Rcpp::export]]
vec SstarGibbs(int N, int M, vec z, vec totalmu, mat sigma){
    //vec zSubvec = z.subvec(0,M-1).t();
    //vec zSub;
    // temp = mvrnormArma(1, z , sigma);
    //mat temp; //= z.subvec(0, M-1)*(z.subvec(0, M-1).t());
    vec mustar(N);
    for(int k=1; k<=(N/M); k++){
      //zSub = z.subvec((k-1)*M, k*M-1);
    	mustar.subvec((k-1)*M, k*M-1) = mvrnormArma(1, sigma*(z.subvec((k-1)*M, k*M-1)-
  							  totalmu.subvec((k-1)*M, k*M-1 ) ), sigma).t(); // mustar[] should be Mx1 vector 
	};
return(mustar) ; }

// [[Rcpp::export]]
mat testfcn(vec vec1){
  mat mat1 = eye(3,3);
  mat mat2 = eye(5,5);
  mat mat3 = kron(mat1, mat2);
  return vec1*mat3*vec1.t();
  } 
  
//******************************************************************
//******************* main function *************************************
// [[Rcpp::export]]
List ZIP_mcmc (
          const double k,
          const int iters,
          const NumericVector y = 0,
          const NumericVector idx_m = 0,
          const NumericVector idx_mt = 0,
          const NumericVector idx_mta = 0,
          const NumericVector idx_mtab = 0,
          const NumericVector idx_mtabc = 0
          ){
    //idx_mtabc = ((c-1)*(M*T*A*B) + (b-1)*(M*T*A) + (a-1)*(M*T) + (t-1)*M + m  ) 
    // if c, b, a, t, m>=1 
    //   c = idx_mtabc / M*T*A*B + 1
    //   b = (idx_mtabc % M*T*A*B ) / (M*T*A) + 1      
    //   a = (idx_mtabc % M*T*A) / (M*T) + 1
    //   t = (idx_mtabc % M*T) / M + 1
    //   m = (idx_mtabc) % M
    //
    // if c, b, a, t, m >=0
    // idx_mtabc = (c*(M*T*A*B) + b*(M*T*A) + a*(M*T) + t*M + m  ) 
    //   c = idx_mtabc / M*T*A*B 
    //   b = (idx_mtabc % M*T*A*B ) / (M*T*A)
    //   a = (idx_mtabc % M*T*A) / (M*T) 
    //   t = (idx_mtabc % M*T) / M 
    //   m = (idx_mtabc) % M
    
    //define inverse wishart
      Environment env("package:MCMCpack");
      Function riwish = env["riwish"];
      
    //constants
        const int M = 24; const int T = 11; const int A = 20; const int B = 7; const int C  = 20; const int N = 2057928;
    
    //priors 
        const int tau_prior_a = 1.0; const int tau_prior_b = 1.0; 
        const int nu_STM = 2*M; const int nu_SAM = 2*M; const int nu_SBM = 2*M; const int nu_SCM = 2*M; 
        const int nu_RC = 2*M;
        const mat Psi_RC = eye(M,M); 
        const mat Psi_STM = eye(M,M); const  mat Psi_SAM = eye(M,M); const mat Psi_SBM = eye(M,M);
        const mat Psi_SCM = eye(M,M);
        mat sigmasqSM = eye(M,M);
        
        
        
    //declare matrices
        int iter = 0;
        int dof = 0; 
        mat Stilde = zeros(M,T); mat lmbda = zeros(M,T);
        int m = 0; int t = 0; int a = 0; int b = 0; int c = 0; 
        int mt = 0; int mta = 0; int mtab = 0; int mtabc =0; 
        int i; //counter  
        vec z(N); vec p(N); vec beta0 = rnorm(M, 0, 1e-4); 
        vec mu = rnorm(N, 0, 1e-5);
        mat RC = eye(M,M); mat RCinv = eye(M,M);
        vec SstarM(N); 
        vec SM = rnorm(M, 0, 1e-3);  
        vec ST = rnorm(M*T, 0, 1e-3);
        vec SA = rnorm(M*T*A, 0, 1e-3);
        vec SB = rnorm(M*T*A*B, 0, 1e-3);
        vec SC = rnorm(M*T*A*B*C, 0, 1e-3);
        vec SstarM_new(N);
        vec SM_new = rnorm(M, 0, 1e-3);
        vec ST_new = rnorm(M*T, 0, 1e-3);
        vec SA_new = rnorm(M*T*A, 0, 1e-3);
        vec SB_new = rnorm(M*T*A*B, 0, 1e-3);
        vec SC_new = rnorm(M*T*A*B*C, 0, 1e-3);
        
        mat STM = eye(M,M); mat SAM = eye(M,M); mat SBM = eye(M,M); mat SCM = eye(M,M);
        mat TTinv = eye(T,T); mat TAinv = eye(T,T); mat TBinv = eye(T,T); mat TCinv = eye(T,T);
        //double temp_prob;
        double sigma_Msq = 1; 
       // double temp_logmu_nostar = 0; double temp_logmu_withstar = 0; 
        double temp = 0;
        mat tempmat = eye(M,M);
        double loglikeli_diff = 0; //double loglikeli_new = 0; double loglikeli_old = 0;
        vec temparr;
        vec temparr2;
        SM = mvrnormArma(1, beta0, sigmasqSM).t(); 
     
       
for (iter=0; iter<iters; iter++){

        //mat TT = generateTmat(0.5);
    //sample p
        for ( i = 1; i<= N; i++){
          temp = pow(mu[i], -k)*exp(k*SstarM[i]);
          p[i] = temp / (1.0 + temp);
        }    
    //sample z  
        //P(z=1 | y=0) = p / (p + (1-p)*exp(-mu))
        //P(z=1 | y=1) = 1
        temparr = runif(N);
        for ( i = 1; i<=N; i++){
           if (y[i] > 0){z[i]=1;} 
           else if( temparr[i] < p[i]/(p[i]+(1-p[i])*(exp(-mu[i])))){z[i] = 1;}
           else {z[i]=0; }      
        }
      
    // RC 
      dof =  nu_RC + (N/M) ;
      tempmat = zeros(M,M);
      for(i=0; i<(N/M); i++){
        temparr = SstarM.subvec(i*M, (i+1)*M-1);
      	tempmat += temparr*(temparr.t()) ; // mustar[] should be Mx1 vector 
	    };
      lmbda = Psi_RC + tempmat;
      RC = as<arma::mat>(riwish(dof, lmbda));
      RCinv = inv(RC);
  
 
    // beta0
        temp = 1.0/(sigma_Msq +1.0);
        beta0 = mvrnormArma(1, SM / temp, 1.0/temp *eye(M,M)).t();
        //beta0 = mvrnormArma(1, SM , eye(M,M)).t();
        
    // sigma_M^2
        temp = R::rgamma(tau_prior_a + 0.5*M, tau_prior_b + 0.5*as_scalar((beta0-SM).t()*(beta0-SM)) );
        sigma_Msq = 1.0/ temp; 
      
    // S^*
         SstarM_new =  SstarM + 1e-2*armaNormal(N);
        loglikeli_diff = 0;
        for (i = 1; i<= N; i ++){
          m = idx_m[i];
          temp = log(mu[i]) - SstarM[i];

          if (y[i] >0){ 
            loglikeli_diff += y[i]*(SstarM_new[i] - SstarM[i])+
                              exp(temp)*(-exp(SstarM_new[i]) + exp(SstarM[i]) );}
          else if (z[i]==0){
                loglikeli_diff += exp(temp)*(-exp(SstarM_new[i]) + exp(SstarM[i]) ); 
                      }
        }
      //update log likelihood difference from the hierarchical part
          for(i=0; i<(N/M); i++){
                temparr = SstarM_new.subvec(i*M, (i+1)*M-1);
                temparr2 = temparr.t()*RCinv*temparr;
                loglikeli_diff += -0.5*temparr2[0] ; 
                
                temparr = SstarM.subvec(i*M, (i+1)*M-1);
                temparr2 = temparr.t()*RCinv*temparr;
                loglikeli_diff += 0.5*temparr2[0] ;  ; 
        	 };
  
        
        temparr = runif(1);
        if (temparr[1] < std::min(1.0, exp(loglikeli_diff))  ){
          SM = SM_new;
          for (i=1; i<=N; i++){mu[i] *= exp(SM_new[m]-SM[m]);} 
        }  

    // S^M
        SM_new =  mvrnormArma(1, SM, 1e-3*eye(M,M)).t();
        loglikeli_diff = 0;
        for (i = 1; i<= N; i ++){
          m = idx_m[i];
          temp = log(mu[i]) - SM[m] - SstarM[i];

          if (y[i] >0){ 
            loglikeli_diff += -log(1+exp(-k*(temp+SM_new[m]))) + 
                               log(1+exp(-k*(temp+SM[m])))+
                                y[i]*(SM_new[m] - SM[m])+
                              exp(temp + SstarM[i])*(-exp(SM_new[m]) + exp(SM[m]) );}
          else {
            if (z[i]==1){
                loglikeli_diff +=  -k*(SM_new[m]-SM[m]) -
                                   log(1+exp(-k*(temp+SM_new[m]))) + 
                                   log(1+exp(-k*(temp+SM[m]))); 
                                   }
            else{
                loglikeli_diff +=  -log(1+exp(-k*(temp+SM_new[m]))) +
                                    log(1+exp(-k*(temp+SM[m])))-
                                    exp(temp + SstarM[i])*(-exp(SM_new[m]) + exp(SM[m]) );
            }
          }
        }
        temparr = runif(1);
        if (temparr[1] < std::min(1.0, exp(loglikeli_diff))  ){
          SM = SM_new;
          for (i=1; i<=N; i++){mu[i] *= exp(SM_new[m]-SM[m]);} 
        }

    // S^T
        ST_new =  ST + 1e-2*armaNormal(M*T);
        loglikeli_diff = 0;
        for (i = 1; i<= N; i ++){
          m = idx_m[i];
          mt = idx_mt[i];
          temp = log(mu[i]) - ST[mt] - SstarM[i];

          if (y[i] >0){ 
            loglikeli_diff += -log(1+exp(-k*(temp + ST_new[mt]))) + 
                               log(1+exp(-k*(temp + ST[mt])))+
                                y[i]*(ST_new[mt] - ST[mt])+
                                exp(temp + SstarM[i])*(-exp(ST_new[mt]) + exp(ST[mt]) );}
          else {
            if (z[i]==1){
                loglikeli_diff +=  -k*(ST_new[mt]-ST[mt]) -
                                   log(1+exp(-k*(temp+ST_new[mt]))) + 
                                   log(1+exp(-k*(temp+ST[mt]))); 
                                   }
            else{
                loglikeli_diff +=  -log(1+exp(-k*(temp+ST_new[mt]))) +
                                    log(1+exp(-k*(temp+ST[mt])))+
                                    exp(temp + SstarM[i])*(-exp(ST_new[mt]) + exp(ST[mt]) );
            }
          }
        }
        //update log likelihood difference from the hierarchical part
            tempmat = ST_new.t()*(kron(inv(STM), TTinv))*ST_new;
            loglikeli_diff += -1/2* tempmat[0]; 
            tempmat = ST.t()*kron(inv(STM), TTinv)*ST;
            loglikeli_diff +=  1/2* tempmat[0];
        
        temparr = runif(1);
        if (temparr[1] < std::min(1.0, exp(loglikeli_diff))  ){
          ST = ST_new;
          for (i=1; i<=N; i++){mu[i] *= exp(ST_new[mt]-ST[mt]);} 
        }
    
    // S^A
        SA_new =  SA + 1e-2*armaNormal(M*T*A);
        loglikeli_diff = 0;
        for (i = 1; i<= N; i ++){
          m = idx_m[i];
          mta = idx_mta[i]; 
          temp = log(mu[i]) - SA[mta] - SstarM[i];

          if (y[i] >0){ 
            loglikeli_diff += -log(1+exp(-k*(temp + SA_new[mta]))) + 
                               log(1+exp(-k*(temp + SA[mta])))+
                                y[i]*(SA_new[mta] - SA[mta])+
                                exp(temp + SstarM[i])*(-exp(SA_new[mta]) + exp(SA[mta]) );}
          else {
            if (z[i]==1){
                loglikeli_diff +=  -k*(SA_new[mta]-SA[mta]) -
                                   log(1+exp(-k*(temp+SA_new[mta]))) + 
                                   log(1+exp(-k*(temp+SA[mta]))); 
                                   }
            else{
                loglikeli_diff +=  -log(1+exp(-k*(temp+SA_new[mta]))) +
                                    log(1+exp(-k*(temp+SA[mta])))-
                                    exp(temp + SstarM[i])*(-exp(SA_new[mta]) + exp(SA[mta]) );
            }
          }
        }
       //update log likelihood difference from the hierarchical part
          for(a=0; a<A; a++){   
            temparr = SA_new.subvec(a*M*T, (a+1)*M*T-1);
            tempmat = temparr.t()*(kron(inv(SAM), TAinv))*temparr;
            loglikeli_diff += -1/2* tempmat[0]; 
            temparr = SA.subvec(a*M*T, (a+1)*M*T-1);
            tempmat = temparr.t()*kron(inv(SAM), TAinv)*temparr;
            loglikeli_diff +=  1/2* tempmat[0];
          }        
        temparr = runif(1);
        if (temparr[1] < std::min(1.0, exp(loglikeli_diff))  ){
          SA = SA_new;
          for (i=1; i<=N; i++){mu[i] *= exp(SA_new[mta]-SA[mta]);} 
        }
        
    // S^B
        SB_new =  SB + 1e-2*armaNormal(M*T*A*B);
        loglikeli_diff = 0;
        for (i = 1; i<= N; i ++){
          m = idx_m[i];
          mtab = idx_mtab[i]; 
          temp = log(mu[i]) - SB[mtab] - SstarM[i];

          if (y[i] >0){ 
            loglikeli_diff += -log(1+exp(-k*(temp+SB_new[mtab]))) + 
                               log(1+exp(-k*(temp+SB[mtab])))+
                                y[i]*(SB_new[mtab] - SB[mtab])+
                                exp(temp + SstarM[i])*(-exp(SB_new[mtab]) + exp(SB[mtab]) );}
          else {
            if (z[i]==1){
                loglikeli_diff +=  -k*(SB_new[mtab]-SB[mtab]) -
                                   log(1+exp(-k*(temp+SB_new[mtab]))) + 
                                   log(1+exp(-k*(temp+SB[mtab]))); 
                                   }
            else{
                loglikeli_diff +=  -log(1+exp(-k*(temp+SB_new[mtab]))) +
                                    log(1+exp(-k*(temp+SB[mtab])))+
                                    exp(temp + SstarM[i])*(-exp(SB_new[mtab]) + exp(SB[mtab]) );
            }
          }
        }
         //update log likelihood difference from the hierarchical part
        for(b=0; b<B;b++){  
            for(a=0; a<A; a++){   
              temparr = SB_new.subvec(b*M*T*A + a*M*T, b*M*T*A +(a+1)*M*T-1);
              tempmat = temparr.t()*(kron(inv(SBM), TBinv))*temparr;
              loglikeli_diff += -1/2* tempmat[0]; 
              temparr = SB.subvec(b*M*T*A + a*M*T, b*M*T*A +(a+1)*M*T-1);
              tempmat = temparr.t()*kron(inv(SBM), TBinv)*temparr;
              loglikeli_diff +=  1/2* tempmat[0];
            }  
        }  
        temparr = runif(1);
        if (temparr[1] < std::min(1.0, exp(loglikeli_diff))  ){
          SB = SB_new;
          for (i=1; i<=N; i++){mu[i] *= exp(SB_new[mtab]-SB[mtab]);} 
        }
  // S^C
        SC_new =  SC + 1e-2*armaNormal(M*T*A*B*C);
        loglikeli_diff = 0;
        for (i = 1; i<= N; i ++){
          m = idx_m[i];
          mtabc = idx_mtabc[i]; 
          temp = log(mu[i]) - SC[mtabc] - SstarM[i];

          if (y[i] >0){ 
            loglikeli_diff += -log(1+exp(-k*(temp+SC_new[mtabc]))) + 
                               log(1+exp(-k*(temp+SC[mtabc])))+
                                y[i]*(SC_new[mtabc] - SC[mtabc])+
                                exp(temp + SstarM[i])*(-exp(SC_new[mtabc]) + exp(SC[mtabc]) );}
          else {
            if (z[i]==1){
                loglikeli_diff +=  -k*(SC_new[mtabc]-SC[mtabc]) -
                                   log(1+exp(-k*(temp+SC_new[mtabc]))) + 
                                   log(1+exp(-k*(temp+SC[mtabc]))); 
                                   }
            else{
                loglikeli_diff +=  -log(1+exp(-k*(temp+SC_new[mtabc]))) +
                                    log(1+exp(-k*(temp+SC[mtabc])))-
                                   exp(temp + SstarM[i])*(-exp(SC_new[mtabc]) + exp(SC[mtabc]) );
            }
          }
        }
        //update log likelihood difference from the hierarchical part
        for(c=0; c<C; c++){  
          for(b=0; b<B;b++){  
              for(a=0; a<A; a++){   
                temparr = SC_new.subvec(c*M*T*A*B + b*M*T*A + a*M*T, 
                                        c*M*T*A*B + b*M*T*A + (a+1)*M*T-1);
                tempmat = temparr.t()*(kron(inv(SCM), TCinv))*temparr;
                loglikeli_diff += -1/2* tempmat[0]; 
                temparr = SC.subvec(c*M*T*A*B + b*M*T*A + a*M*T, 
                                        c*M*T*A*B + b*M*T*A + (a+1)*M*T-1);
                tempmat = temparr.t()*kron(inv(SCM), TCinv)*temparr;
                loglikeli_diff +=  1/2* tempmat[0];
              }  
            }  
          }
          temparr = runif(1);
         if (temparr[1] < std::min(1.0, exp(loglikeli_diff))  ){
          SC = SC_new;
          for (i=1; i<=N; i++){mu[i] *= exp(SC_new[mtabc]-SC[mtabc]);} 
        }

    
  //  S^T_M
      dof = T + nu_STM;
      temparr = ST;
      lmbda = Psi_STM;
      Stilde = zeros(M,T);
      for(t=0; t<T; t++){
         Stilde.col(t) = temparr.subvec(t*M, (t+1)*M-1);
      }
      lmbda += Stilde * TTinv * (Stilde.t()); 
      STM = as<arma::mat>(riwish(dof, lmbda));
  
  //S^A_M
      temparr = SA;
      dof = T*A +nu_SAM;
      lmbda = Psi_SAM;
      for (a =0; a<A; a++){
         //define Stilde
         for(t=0; t<T; t++){
         Stilde.col(t) = temparr.subvec(a*M*T + t*M , a*M*T + (t+1)*M-1);
        }
        lmbda +=  Stilde * TAinv * (Stilde.t()); 
      }
      SAM = as<arma::mat>(riwish(dof, lmbda));
  
  //S^B_M
      temparr = SB;
      dof = T*A*B +nu_SBM;
      lmbda = Psi_SBM;
      for (b =0; b<B; b++){
        for (a =0; a<A; a++){
           //define Stilde
           for(t=0; t<T; t++){
           Stilde.col(t) = temparr.subvec(b*M*T*A + a*M*T + t*M , 
                                          b*M*T*A + a*M*T + (t+1)*M-1);}
           lmbda +=  Stilde * TAinv * (Stilde.t()); 
        }
      }
      SBM = as<arma::mat>(riwish(dof, lmbda));
  //S^C_M
      temparr = SC;
      dof = T*A*B*C +nu_SCM;
      lmbda = Psi_SCM;
      for (c =0; c<C; c++){
        for (b =0; b<B; b++){
          for (a =0; a<A; a++){
             //define Stilde
             for(t=0; t<T; t++){
             Stilde.col(t) = temparr.subvec(c*M*T*A*B + b*M*T*A + a*M*T + t*M , 
                                            c*M*T*A*B + b*M*T*A + a*M*T + (t+1)*M-1);}
             lmbda +=  Stilde * TAinv * (Stilde.t()); 
          }
        }
      }
      SCM = as<arma::mat>(riwish(dof, lmbda));
  
  if (iter%100==0){
      Rcpp::Rcout << "iteration = " << iter << std::endl;
    }
  } //end of iterations 
        
        
// returns
    List ret ;
   // ret["TT"] = TT ;
    ret["RC"] = RC ; 
    ret["SM"] = SM ; 
    ret["STM"] = STM ; 
    ret["SAM"] = SAM ; 
    ret["SBM"] = SBM ;
    ret["SCM"] = SCM ;
    //ret["lmbda"] = lmbda ; 
    ret["beta0"] = beta0 ;
    //ret["prob"] = p ;
    //ret["z"] = z;
  return(ret) ;
}





/*  
// [[Rcpp::export]]
NumericMatrix inverseWish(int Dof, NumericMatrix lmbda){
  Environment myEnv = Environment::global_env();
  //Environment stats("package:MCMCpack");
  Function riwish = myEnv["riwish"];
  //int dof = Rcpp::as<int>(Dof);
  //NumericMatrix lmbda(Lmbda);
  NumericMatrix Ret = riwish(Named("dof", Dof), Named("lmbda", lmbda));
}

// [[Rcpp::export]]
RObject InvWish(int dof, NumericMatrix lmbda, Function f){
  return f(dof, lmbda);
}

// [[Rcpp::export]]
List InvWish2(int dof, NumericMatrix lmbda, Function f){
  List result;
  result["res"] = f(dof, lmbda);
  return result;
}

// [[Rcpp::export]]
NumericMatrix InvWish3(int dof, NumericMatrix lmbda){
  //List result;
  Environment env = Environment::global_env();
  Function riwish = env["riwish"];
  //result["res"] = riwish(dof, lmbda);
  return  riwish(dof, lmbda);
}


// [[Rcpp::export]]
NumericMatrix InvWish4(int dof, NumericMatrix lmbda){
  //List result;
  Environment env("package:MCMCpack");
  Function riwish = env["riwish"];
  //result["res"] = riwish(dof, lmbda);
  return riwish(dof, lmbda);
}
*/ 

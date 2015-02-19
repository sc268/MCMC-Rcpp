#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat RClambda (int N, int M, NumericVector mustar) {
    vec mustarvec = as<vec>(mustar);
    vec mustarSubvec;
    mat temp = mustarvec.subvec(0, M-1)*(mustarvec.subvec(0, M-1).t());
    for(int k=2; k<=(N/M); k++){
    	mustarSubvec = mustarvec.subvec((k-1)*M, k*M-1);
  		temp += mustarSubvec*(mustarSubvec.t()); // mustar[] should be Mx1 vector 
	};
return(temp) ; }



#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector RC_lambda (int N, int M, NumericVector mustar) {
	NumericVector temp;
	for(int k=0; k<(N/M); k++){
		temp += dot(mustar[(1+(k-1)*M):(k*M)],mustar[(1+(k-1)*M):(k*M)]); 
		// mustar[] should be Mx1 vector 
	};
    return(temp) ;
}


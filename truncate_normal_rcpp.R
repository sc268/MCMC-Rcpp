library(Rcpp)

cppFunction("
NumericVector logabs(NumericVector x) { return log(abs(x));
}
            ")
logabs(seq(-5, 5, by=2))


lower = rep(0, N)
upper = rep(0, N)
for(k in 1:N){
  if (Abundance[k]>0){lower[k]=0;    upper[k]=500}else {
                      lower[k]=-500; upper[k]=0}
}
getwd()

library(Rcpp)
sourceCpp('tn.cpp')
sourceCpp('tnALL.cpp')
#get1TN(5,1,0,1000)



ptm <- proc.time()
getTN(N, rnorm(N,0,1), 1, lower, upper)
proc.time() - ptm    




sourceCpp('piSugar.cpp')
set.seed(42); a <- piR(1.0e7) 
set.seed(42); b <- piSugar(1.0e7)
identical(a,b)

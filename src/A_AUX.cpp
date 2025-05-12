#include "A_AUX.h"

double LogSumExp_cpp(arma::rowvec logX){
double a = max(logX);
return(  a + log(accu( exp( logX-a ) )));
}

arma::colvec reverse_cumsum_cpp(arma::colvec X){
  return( accu(X) - arma::cumsum(X));
}

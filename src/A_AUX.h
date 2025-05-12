#ifndef A_AUX
#define A_AUX

#include <RcppArmadillo.h>

double LogSumExp_cpp(arma::rowvec logX);

arma::colvec reverse_cumsum_cpp(arma::colvec X);

#endif


#ifndef C_NORM_CONST
#define C_NORM_CONST

#include <RcppArmadillo.h>

double lbeta_normconst_cpp(double a,
                           double b);

arma::colvec lbeta_normconst_vec_cpp(arma::colvec a,
                                     arma::colvec b);

arma::mat lbeta_normconst_mat_cpp(arma::mat a, 
                                  arma::mat b);

arma::colvec ldirichlet_normconst_vec_cpp(arma::mat beta_lk);

double ldirichlet_normconst_cpp(arma::colvec alpha_k);


#endif


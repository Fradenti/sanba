#include "C_NORM_CONST.h"

double lbeta_normconst_cpp(double a, double b){
  double C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - (C_ab));
}

arma::colvec lbeta_normconst_vec_cpp(arma::colvec a, arma::colvec b){
  arma::colvec C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - (C_ab));
}

arma::mat lbeta_normconst_mat_cpp(arma::mat a, arma::mat b){
  arma::mat C_ab = lgamma(a) + lgamma(b) - lgamma(a+b);
  return( - C_ab);
}


arma::colvec ldirichlet_normconst_vec_cpp(arma::mat beta_lk){
  
  int K = beta_lk.n_cols;
  arma::colvec C(K);
  for(int k=0; k<K; k++){
    C(k) = lgamma(arma::accu(beta_lk.col(k))) -
                  arma::accu(lgamma(beta_lk.col(k)));
  }
  
  return(C);
}



double ldirichlet_normconst_cpp(arma::colvec alpha_k){
  
  double  C = lgamma(arma::accu(alpha_k)) -
      arma::accu(lgamma(alpha_k));
  
  return(C);
}

#include "B_EXP_VALUES.h"

// slightly faster
arma::colvec E_log_beta(arma::colvec a,
                        arma::colvec b){
  
  int n = a.n_elem;
  arma::colvec res(n);
  for(int i=0; i<n; i++){
    res[i] = R::digamma(a[i]) - R::digamma(a[i] + b[i]);
  }
  return(res);
}

// -----------------------------------------------------------------------------

arma::colvec E_log_IG(arma::colvec a,
                      arma::colvec b){
  // convert from arma to vanilla rcpp so I can use vectorized digamma
  Rcpp::NumericVector av  = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(a));
  Rcpp::NumericVector res = digamma(av);
  
  // back to armadillo
  arma::colvec ares = log(b) - Rcpp::as<arma::colvec>(Rcpp::wrap(res));
  return(ares);
}

// -----------------------------------------------------------------------------


//slightly faster
arma::colvec E_log_DIR(arma::colvec a){
  
  int n = a.n_elem;
  arma::colvec res(n);
  double Sum_a = arma::accu(a);
  
  for(int i = 0; i < n; i++){
   res[i] = R::digamma( a[i] );   
  }
  
  return(res - R::digamma(Sum_a));
}

// -----------------------------------------------------------------------------

arma::colvec E_log_p_Y_Mtheta_cpp(arma::colvec Y,
                                  double ml,
                                  double kl,
                                  double al,
                                  double bl){
  return(
    - .5 * (
      log(bl) - R::digamma(al)   +
        1/kl + ( (al/bl) * (Y-ml) % (Y-ml) )  )
  );
}



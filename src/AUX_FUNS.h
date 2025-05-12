#ifndef AUX_FUNS
#define AUX_FUNS

#include <RcppArmadillo.h>

arma::vec rdirichlet(arma::vec par) ;

int sample_i(arma::vec ids, arma::vec prob) ;

arma::vec relabel_arma(arma::vec cluster) ;

Rcpp::List sample_model_parameters(const arma::vec& y, 
                                   arma::vec M_iter,
                                   int maxL,
                                   double m0, double tau0,
                                   double lambda0, double gamma0) ;

double fun_xi(double kappa, int i) ;

int compute_trunc(double u_min, double kappa) ;

arma::vec stick_breaking(arma::vec beta_var);

#endif


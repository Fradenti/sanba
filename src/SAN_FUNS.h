#ifndef SAN_FUNS
#define SAN_FUNS

#include "AUX_FUNS.h"

arma::vec dirichlet_sample_distr_weights(arma::vec S_iter, 
                                         double alpha, 
                                         int maxK);

arma::mat dirichlet_sample_obs_weights(arma::vec M_iter, 
                                       arma::vec clusterD_long,
                                       double beta,
                                       int K_iter, int L_iter,
                                       int maxK, int maxL);


arma::vec sample_distr_cluster(const arma::vec& y, const arma::vec& group,
                               arma::vec pi, arma::mat omega,
                               int maxK, int maxL,
                               arma::vec mu, arma::vec sigma2);

arma::vec sample_obs_cluster(const arma::vec& y,
                             arma::vec clusterD_long,
                             arma::mat omega,
                             int maxL,
                             arma::vec mu, arma::vec sigma2);

Rcpp::List fisan_weights_update_slice_sampler(const arma::vec& group, 
                                              arma::vec S_iter, 
                                              double alpha,
                                              arma::vec xi,
                                              int maxK);

#endif
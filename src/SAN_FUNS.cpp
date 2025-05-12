#include "SAN_FUNS.h"


/*
 *  UPDATE MIXTURE WEIGHTS
*/
/* sample distributional probabilities pi
 * (pi_1,...,pi_K) | . ~ Dir_K(alpha/K + N_1, ..., alpha/K + N_K)
 * NOTE: we consider a "dynamic mixture", with alpha/K
 * where N_k is the number of groups assigned to distributional cluster k
*/
arma::vec dirichlet_sample_distr_weights(arma::vec S_iter,
                                    double alpha,
                                    int maxK)
{
  // K_iter is the dimension of the Dirichlet rv
  // maxK is the size of the allocated vector, and I need to fill only the first K_iter elements
  arma::vec out_pi(maxK, arma::fill::zeros) ;

  /* sample distributional probabilities pi */
  arma::vec dir_paramD(maxK) ;
  for(int k = 0; k < maxK ; k++)
  {
    arma::uvec ind_k = find(S_iter == k) ;
    dir_paramD(k) = alpha + ind_k.n_elem ;
  }
  out_pi(arma::span(0, maxK-1)) = rdirichlet(dir_paramD) ;

  return(out_pi) ;
}


/*
 * for each distributional cluster k, sample observational probabilities omega
 * (om_1k,...,om_Lk) | . ~ Dir_L(beta/L + n_1k, ..., beta/L + n_Lk)
 * NOTE: we consider a "dynamic mixture", with beta/L
 * where n_lk is the number of observations assigned to cluster l in the distributional cluster k
*/
arma::mat dirichlet_sample_obs_weights(arma::vec M_iter,
                                       arma::vec clusterD_long,
                                       double beta,
                                       int K_iter, int L_iter,
                                       int maxK, int maxL)
{
  // as before,
  // L_iter is the dimension of the observational Dirichlet rv, K_iter the current dimension of the distrib. rv
  // maxL x maxK is the size of the allocated matrix, and I need to fill only the first L_iter x K_iter elements

  arma::mat out_omega(maxL, maxK, arma::fill::zeros) ;

  /* for each distributional cluster k, sample observational probabilities omega */
  arma::vec dir_param(L_iter) ;
  for(int k = 0; k < K_iter; k++)
  {
    out_omega.col(k).fill(0) ;
    arma::uvec ind_clusterD_k = find(clusterD_long == k) ;
    arma::vec subcluster = M_iter.elem( ind_clusterD_k ) ;

    dir_param.fill(0) ;
    for(int l = 0; l < L_iter ; l++)
    {
      arma::uvec subcluster_l = find( subcluster == l ) ;
      dir_param(l) = beta + subcluster_l.n_elem ; // beta/L_iter + subcluster_l.n_elem ;
    }
    out_omega(arma::span(0, L_iter-1), arma::span(k)) = rdirichlet(dir_param) ;
  }

  return(out_omega) ;
}


arma::vec sample_distr_cluster(const arma::vec& y, const arma::vec& group,
                               arma::vec pi, arma::mat omega,
                               int maxK, int maxL,
                               arma::vec mu, arma::vec sigma2)
{
  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ;

  arma::vec distr_cluster_id = arma::linspace(0, maxK-1, maxK) ;
  arma::vec probD(maxK) ;

  arma::vec out(G) ;

  // S_j is categorical, with Pr( S_j = k | - ) = pi_k * ( om_1k^#(M_ij = 1) ... om_Lk^#(M_ij = L) )
  // hence the log is:  log(pi_k) + sum_l  #(M_ij = l) * log( om_lk )
  for(int j = 0; j < G; j++)
  {
    arma::uvec ind_group_j = find(group == j) ;  // restrict to observations in the j-th population
    arma::vec mixdens(ind_group_j.n_elem) ;

    for(int k = 0; k < maxK; k++)  // I have to compute the prob for each k = 1,..., maxK
    {
      for(unsigned int i = 0; i < ind_group_j.n_elem; i++)
      {
        mixdens(i) = 0 ;
        for(int l = 0; l < maxL; l++) {
          mixdens(i) = mixdens(i) + exp( log(omega(l, k)) + R::dnorm( y(ind_group_j(i)), mu(l), std::sqrt(sigma2(l)), true ) );
        }
        mixdens(i) = log( mixdens(i) ) ;
      }
      probD(k) =  log( pi(k) ) + arma::accu(mixdens) ;
    }
    probD =  probD - max(probD) ;
    probD = exp(probD) ;
    out(j) = sample_i(distr_cluster_id, probD) ;
  }

  return(out) ;
}



arma::vec sample_obs_cluster(const arma::vec& y,
                             arma::vec clusterD_long,
                             arma::mat omega,
                             int maxL,
                             arma::vec mu, arma::vec sigma2)
{
  arma::vec obs_cluster_id = arma::linspace(0, maxL-1, maxL) ;
  int N = y.n_elem ;
  arma::vec probO(maxL) ;

  arma::vec out(N) ;

  // M_ij is categorical, with Pr(M_ij = l | S_j = k) = om_lk * Pr(y_ij | mu_l)
  for(int i = 0; i < N; i++)
  {
    for(int l = 0; l < maxL; l++) {
      probO(l) = log( omega( l, clusterD_long(i) ) ) + R::dnorm(y(i), mu(l), std::sqrt(sigma2(l)), true) ;
    }
    probO =  probO - max(probO) ;
    probO = exp(probO) ;
    out(i) = sample_i(obs_cluster_id, probO) ;
  }

  return(out) ;
}

Rcpp::List fisan_weights_update_slice_sampler(const arma::vec& group,
                                              arma::vec S_iter,
                                              double alpha,
                                              arma::vec xi,
                                              int maxK)
{
  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ;
  arma::vec out_pi(maxK, arma::fill::zeros) ;

  int maxK_new ;
  arma::vec u_D(G) ;

  int a_k ; double b_k ;

  /*---------------------------------------------*/
  /*
   *  SAMPLE SLICE - uniform r.v.
   */
  /* distributional */
  for(int j = 0; j < G; j++)
  {
    u_D(j) = R::runif( 0, xi(S_iter(j)) ) ;
  }
  int maxK_slice_tmp = compute_trunc(u_D.min(), 0.5) + 1 ;
  maxK_new = std::min(maxK, maxK_slice_tmp) ;

  /*---------------------------------------------*/


  /*---------------------------------------------*/
  /*
   *  UPDATE MIXTURE WEIGHTS
   */
  /* sample distributional probabilities pi */
  /* (pi_1,...,pi_maxK) | . ~ DP */

  arma::vec v_k(maxK_new) ;
  arma::vec pi_k(maxK_new) ;
  for(int k = 0; k < maxK_new ; k++)
  {
    arma::uvec ind_k = find(S_iter == k) ;
    arma::uvec ind_mk = find(S_iter > k) ;

    a_k = 1 + ind_k.n_elem ;
    b_k = alpha + ind_mk.n_elem ;
    v_k(k) = R::rbeta(a_k, b_k) ;
  }
  pi_k = stick_breaking( v_k ) ;
  out_pi(arma::span(0, maxK_new-1)) = pi_k ;


  /*---------------------------------------------*/

  return Rcpp::List::create(Rcpp::Named("new_maxK") = maxK_new,
                            Rcpp::Named("new_pi") = out_pi,
                            Rcpp::Named("u_D") = u_D);

}

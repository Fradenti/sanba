#include "CAM_FUNS.h"






// This function performs steps 1-4 of Algorithm 1 of Denti et al. (2021):
// - sample observational and distributional slice variables,
// - update mixture weights pi and omega

// The function outputs:
// - the number of needed mixture components (observational and distributional)
// - distributional and observational probabilities (pi and omega)
// - the observational beta r.v. which are then necessary for the update the DP parameter
// - uniform "slice" random variables (necessary for the update of the cluster allocation)
Rcpp::List weights_update_slice_sampler(const arma::vec& y, const arma::vec& group, 
                         arma::vec S_iter, arma::vec M_iter, 
                         arma::vec clusterD_long,
                         arma::vec xi,
                         double & alpha, double & beta, 
                         int & maxK, int & maxL)
{
  int N = y.n_elem ;
  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ;
  
  arma::vec out_pi(maxK, arma::fill::zeros) ;
  arma::mat out_omega(maxL, maxK, arma::fill::zeros) ;
  
  int maxK_new ; int maxL_new ;
  arma::vec u_O(N) ; arma::vec u_D(G) ;
  
  int a_k ; double b_k ;
  int a_lk ; double b_lk ;
  
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
  
  /* observational */
  for(int i = 0; i < N; i++)
  {
    u_O(i) = R::runif( 0, xi(M_iter(i)) ) ;
  }
  int maxL_slice_tmp = compute_trunc(u_O.min(), 0.5) + 1 ;
  maxL_new = std::min(maxL, maxL_slice_tmp) ;
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
  
  
  /* for each distributional cluster k, sample observational probabilities omega */
  /* (om_1k,...,om_Lk) | . ~ DP(beta + n_1k, ..., beta + n_Lk)
   where n_lk is the number of observations assigned to cluster l in the distributional cluster k */
  arma::vec v_lk(maxL_new) ;
  arma::vec om_lk(maxL_new) ;
  arma::mat obs_beta_rv(maxL_new, maxK_new) ;
  
  for(int k = 0; k < maxK_new ; k++)
  {
    arma::uvec ind_clusterD_k = find(clusterD_long == k) ; // indices of the y_ij s.t. S_j = k
    arma::vec subcluster = M_iter( ind_clusterD_k ) ; // obs cluster label M_ij of the y_ij ~ G*_k
    
    for(int l = 0; l < maxL_new; l++)
    {
      a_lk = 1;
      b_lk = beta;
      
      if( subcluster.n_elem > 0 )
      {
        arma::uvec ind_lk = find(subcluster == l) ;
        arma::uvec ind_mlk = find(subcluster > l) ;
        
        a_lk = 1 + ind_lk.n_elem ;
        b_lk = beta + ind_mlk.n_elem ;
      }
      v_lk(l) = R::rbeta(a_lk, b_lk) ;
    }
    obs_beta_rv.col(k) = v_lk ;
    om_lk = stick_breaking( v_lk ) ;
    out_omega(arma::span(0, maxL_new-1), arma::span(k)) = om_lk ;
    
  }
  /*---------------------------------------------*/
  
  return Rcpp::List::create(Rcpp::Named("new_maxK") = maxK_new,
                            Rcpp::Named("new_maxL") = maxL_new,
                            Rcpp::Named("new_pi") = out_pi,
                            Rcpp::Named("new_omega") = out_omega,
                            Rcpp::Named("obs_beta_rv") = obs_beta_rv,
                            Rcpp::Named("u_O") = u_O,
                            Rcpp::Named("u_D") = u_D);
  
}




// This function performs step 5 of Algorithm 1 of Denti et al. (2021):
// sample distributional cluster allocation S_j
arma::vec slicedDP_sample_distr_cluster(const arma::vec& group,
                                        arma::vec M_iter,
                                        arma::vec pi, arma::mat omega,
                                        arma::vec u_D,
                                        arma::vec xi,
                                        int maxK_iter)
{

  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ;
  
  arma::vec distr_cluster_id = arma::linspace(0, maxK_iter-1, maxK_iter) ;
  arma::vec probD(maxK_iter) ;
  arma::vec out(G) ;
  double sumdens ;

  // S_j is categorical, with Pr( S_j = k | - ) = I(u_j < xi_k) pi_k/xi_k * ( om_1k^#(M_ij = 1) ... om_Lk^#(M_ij = L) )
  // hence the log is:  log(pi_k) - log(xi_k) + sum_l  #(M_ij = l) * log( om_lk )
  for(int j = 0; j < G; j++)
  {
    arma::uvec ind_group_j = find(group == j) ;  // restrict to observations in the j-th population
    arma::vec mixdens(ind_group_j.n_elem) ;
   
    for(int k = 0; k < maxK_iter; k++)  // I have to compute the prob for each k = 1,..., K_iter
    {
      for(unsigned int i = 0; i < ind_group_j.n_elem; i++) { mixdens(i) = log( omega(M_iter(ind_group_j(i)), k) ) ; }
      sumdens = arma::accu(mixdens) ;
      if(!arma::is_finite(sumdens)) { sumdens = log(0) ;}
      probD(k) =  log( pi(k) ) - log( xi(k) ) + sumdens + log(u_D(j) < xi(k));
    }
    double tmp = max(probD) ;
    for(int k = 0; k < maxK_iter; k++) {  probD(k) = probD(k) - tmp ; }
    probD = exp(probD) ;
    out(j) = sample_i(distr_cluster_id, probD) ;
  }
  
  return(out) ;
}

// This function performs step 6 of Algorithm 1 of Denti et al. (2021):
// sample observational cluster allocation M_ij
arma::vec slicedDP_sample_obs_cluster(const arma::vec& y, 
                                      arma::vec clusterD_long,
                                      arma::mat omega,
                                      arma::vec u_O,
                                      arma::vec xi,
                                      int maxL_iter,
                                      arma::vec mu, arma::vec sigma2)
{
  arma::vec obs_cluster_id = arma::linspace(0, maxL_iter-1, maxL_iter) ;
  int N = y.n_elem ;
  arma::vec probO(maxL_iter) ;
  
  arma::vec out(N) ;
  
  // M_ij is categorical, with Pr(M_ij = l | S_j = k) = I(u_i < xi_l) * om_lk / xi_l * Pr(y_ij | mu_l)
  // the logarithm hence is   log( om_lk ) - log( xi_l ) + log( Pr(y_ij | mu_l) )
  for(int i = 0; i < N; i++) 
  {
    for(int l = 0; l < maxL_iter; l++) {
      probO(l) = log(u_O(i) < xi(l)) + 
        log( omega( l, clusterD_long(i) ) ) - log( xi(l) ) + R::dnorm(y(i), mu(l), std::sqrt(sigma2(l)), true) ;
    }
    probO =  probO - max(probO) ;
    probO = exp(probO) ;
    out(i) = sample_i(obs_cluster_id, probO) ;
  }
  
  return(out) ;
}




// The following functions perform step 8 of Algorithm 1 of Denti et al. (2021):
// update distributional DP parameter 
double sample_alpha(double old_alpha,
                    double hyp_alpha1, double hyp_alpha2,
                    int G , arma::vec S_iter)
{
  double out ;
  arma::vec uniques = arma::unique(S_iter) ;
  int distinct_k = uniques.n_elem;
  double log_eta = log( R::rbeta(old_alpha + 1, G) ) ;
  
  double Q = (hyp_alpha1 + distinct_k - 1) / (G * (hyp_alpha2 - log_eta)) ;
  double pi_eta = Q / (1 + Q) ;
  double u = R::runif(0,1) ;
  
  if (u < pi_eta) {
    out  = R::rgamma( hyp_alpha1 + distinct_k, 1/(hyp_alpha2 - log_eta) ) ;
  } else {
    out = R::rgamma( hyp_alpha1 + distinct_k - 1, 1/(hyp_alpha2 - log_eta) ) ;
  }
  
  return(out) ;
}


// update observational DP parameter
double sample_beta(double hyp_beta1, double hyp_beta2,
                   arma::vec S_iter, arma::vec M_iter, 
                   arma::vec clusterD_long,
                   arma::mat obs_beta_rv)
{
  double out ;
  int bar_S = S_iter.max() ;
  int bar_M ;
  double sum_log_mbeta = 0 ;
  
  arma::vec bar_Ms(bar_S, arma::fill::zeros) ;
  
  for(int k = 0; k < bar_S; k++) 
  { 
    arma::uvec ind_clusterD_k = find(clusterD_long == k) ; 
    arma::vec subM = M_iter(ind_clusterD_k) ;
    if( subM.n_elem > 0 ) 
    { 
      bar_Ms(k) = subM.max() ; 
      for(int l = 0; l < subM.max(); l++) 
      {
        sum_log_mbeta += log( 1-obs_beta_rv(l,k) ) ; 
      }
    }
  }
  bar_M = arma::accu(bar_Ms) ;
  
  out = R::rgamma(hyp_beta1 + bar_M, 1/(hyp_beta2 - sum_log_mbeta) ) ;
  return(out) ;
}

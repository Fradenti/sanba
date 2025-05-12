#include "CAM_FUNS.h"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
Rcpp::List sample_CAM_cpp(int nrep, // number of replications of the Gibbs sampler
                           int burn,
                           const arma::vec & y, // input data
                           const arma::vec & group, // group assignment for each observation in the vector y
                           int maxK, // maximum number of distributional clusters - for memory allocation
                           int maxL, // maximum number of observational clusters - for memory allocation
                           double m0, double tau0, // hyperparameters on the N-iG prior on the mean parameter, mu|sigma2 ~ N(m0, sigma2 / tau0)
                           double lambda0, double gamma0, // hyperparameters on the N-iG prior on the variance parameter, 1/sigma2 ~ Gamma(lambda0, gamma0)
                           bool fixed_alpha, bool fixed_beta, // do you want fixed alpha or beta?
                           double alpha, double beta, // if alpha or/and beta fixed
                           double hyp_alpha1, double hyp_alpha2, // hyperparameters of alpha ( par of the distributional DP ) - if alpha not fixed
                           double hyp_beta1, double hyp_beta2,// hyperparameters of beta ( par of the observational DP ) - if beta not fixed
                           arma::vec mu_start, 
                           arma::vec sigma2_start,
                           arma::vec M_start,
                           arma::vec S_start,
                           double alpha_start,
                           double beta_start,
                           bool progressbar ) 
{
  int N = y.n_elem ;
  arma::vec unique_groups = arma::unique(group) ;
  int G = unique_groups.n_elem ; 
  
  //Rcpp::Rcout << nrep-burn;
  
  // allocate output matrices
  arma::mat mu(maxL, nrep - burn, arma::fill::zeros) ; 
  arma::mat sigma2(maxL, nrep - burn, arma::fill::ones) ; 
  arma::mat out_S(G, nrep - burn) ; // distributional clusters
  arma::mat out_M(N, nrep - burn) ; // observational clusters
  arma::mat out_pi(maxK, nrep - burn, arma::fill::zeros) ;
  arma::cube out_omega(maxL, maxK, nrep - burn, arma::fill::zeros) ; 
  arma::vec out_alpha(nrep - burn) ; // DP parameter for the distributional weights 
  arma::vec out_beta(nrep - burn) ; // DP parameter for the observational weights 
  arma::vec out_maxK(nrep - burn) ; 
  arma::vec out_maxL(nrep - burn) ;
  
  // initialization
  arma::vec current_mu = mu_start ;
  arma::vec current_sigma2 = sigma2_start ;
  arma::vec current_M = M_start;
  arma::vec current_S = S_start;
  arma::vec current_pi(maxK) ; 
  current_pi.fill(1.0/(maxK));
  arma::mat current_omega(maxL,maxK) ; current_omega.fill(1.0/(maxL));
  int current_maxK = maxK ;
  int current_maxL = maxL ;
  /*
    mu.col(0) = mu_start ;
  sigma2.col(0) = sigma2_start ;
  out_M.col(0) = M_start ;
  out_S.col(0) = S_start ;
  pi.col(0).fill(1/(maxK*1.0)) ;
  omega.slice(0).fill(1/(maxL*1.0)) ;
  */
    
  double current_beta = 0;
  double current_alpha = 0;
  
  if(fixed_alpha) { 
    out_alpha.fill(alpha) ; 
    current_alpha = alpha_start ;
  } else {
    current_alpha = alpha_start ;
  }
  
  if(fixed_beta) { 
    out_beta.fill(beta) ; 
    current_beta = beta_start;
  } else {
    current_beta = beta_start ;
  }
  
  
  
  // auxiliary quantities
  arma::vec clusterD_long(N) ; 
  for(int i = 0; i < N; i++) { clusterD_long(i) = S_start( group(i) ) ; }
  Rcpp::List out_params ;
  
  // slice variables
  Rcpp::List weights_slice_sampler ;
  
  arma::vec xi( std::max(maxL, maxK) ) ;
  for(int l = 0; l < std::max(maxL, maxK); l++) { xi(l) = fun_xi(0.5, l+1) ; }
  
  bool warningL = false ; arma::vec iters_maxL(nrep) ; int countL = 0 ;
  bool warningK = false ; arma::vec iters_maxK(nrep) ; int countK = 0 ;
  
  // progress bar 
  bool display_progress = progressbar ;
  Progress p(nrep, display_progress) ;
  
  
  /* START */
  for(int iter = 0; iter < nrep ; iter++)
  {
    if( Progress::check_abort() ) { return -1.0 ; }
    p.increment();
      
      
    /*---------------------------------------------*/
    /*
    *  UPDATE MIXTURE WEIGHTS
    */
        
    weights_slice_sampler = weights_update_slice_sampler(y, group,
                                                         current_S, current_M,
                                                         clusterD_long,
                                                         xi,
                                                         current_alpha, 
                                                         current_beta,
                                                         maxK, maxL) ;
      
    arma::vec tmp_pi = weights_slice_sampler["new_pi"] ;
    current_pi = tmp_pi ;
      
    arma::mat tmp_omega = weights_slice_sampler["new_omega"] ;
    current_omega = tmp_omega ;
      
    int maxK_new = weights_slice_sampler["new_maxK"] ;
    current_maxK = maxK_new ;
    if(current_maxK >= maxK) { 
      warningK = true ;
      iters_maxK(countK) = (iter + 1) ;
      countK++ ;
      current_maxK = maxK-1 ; }
      
    int maxL_new = weights_slice_sampler["new_maxL"] ;
    current_maxL = maxL_new ;
    if(current_maxL >= maxL) { 
      warningL = true ;
      iters_maxL(countL) = (iter + 1) ;
      countL++ ;
      current_maxL = maxL-1 ; }
      
    arma::mat obs_beta_rv = weights_slice_sampler["obs_beta_rv"] ;
    arma::vec u_O = weights_slice_sampler["u_O"] ;
    arma::vec u_D = weights_slice_sampler["u_D"] ;
      
    /*---------------------------------------------*/
        
    /*---------------------------------------------*/
    /*
    *  UPDATE CLUSTER ASSIGNMENT
    */
        
    /* update distributional clusters S */
    current_S = slicedDP_sample_distr_cluster(group,
                                                current_M,
                                                current_pi, 
                                                current_omega,
                                                u_D, xi,
                                                current_maxK) ;
    //out_S.col(iter+1) = out_S.col(iter+1) ;
    for(int i = 0; i < N; i++) { clusterD_long(i) = current_S( group(i)) ; }
    
    /* update observational clusters M */
    current_M = slicedDP_sample_obs_cluster(y,
                                            clusterD_long,
                                            current_omega,
                                            u_O, xi,
                                            current_maxL,
                                            current_mu, 
                                            current_sigma2) ;
    
    /*---------------------------------------------*/
    
      
    /*---------------------------------------------*/
    /*
    *  UPDATE PARAMETERS
    */
    out_params = sample_model_parameters(y, current_M,
                                          maxL,
                                          m0, tau0, lambda0, gamma0) ;
    
    arma::vec tmp_mu = out_params["out_mu"] ;
    current_mu = tmp_mu ; 
    
    arma::vec tmp_sigma2 = out_params["out_sigma2"] ;
    current_sigma2 = tmp_sigma2 ; 
    
    /*---------------------------------------------*/
      
      
    /*---------------------------------------------*/
    /*
    *  UPDATE DIRICHLET HYPER-PARAMETERS 
    */
    
    /* sample alpha */
      if(!fixed_alpha) {
        current_alpha = sample_alpha(current_alpha,
                                     hyp_alpha1, hyp_alpha2,
                                     G, current_S) ;
      }
    
    /* sample beta */
    if(!fixed_beta) {
        current_beta = sample_beta(hyp_beta1, hyp_beta2,
                                   current_S, 
                                   current_M,
                                   clusterD_long,
                                   obs_beta_rv) ;
    }
    
    /*---------------------------------------------*/
      
    // saving
        
    if(iter >= burn){
      int ind = iter-burn;
      //Rcpp::Rcout << ind << "\n";
        
      mu.col(ind) = current_mu;
      sigma2.col(ind) = current_sigma2;
      //Rcpp::Rcout << "a";
      out_maxL(ind) = current_maxL;
      out_maxK(ind) = current_maxK;        
      //Rcpp::Rcout << "b";
        
      out_alpha(ind) = current_alpha;  
      out_beta(ind)  = current_beta;  
      //Rcpp::Rcout << "c";
        
      out_pi.col(ind) = current_pi;
      out_omega.slice(ind) = current_omega;
      //      Rcpp::Rcout << "d";
        
      out_S.col(ind) = current_S;
      out_M.col(ind) = current_M;
  //      Rcpp::Rcout << "e";
        
     }
   }
  
  Rcpp::List warnings ;
  if(warningK & !(warningL)) {
    warnings = Rcpp::List::create(Rcpp::Named("top_maxK") = iters_maxK.head(countK) );
  }
  if(warningL & !(warningK)) {
    warnings = Rcpp::List::create(Rcpp::Named("top_maxL") = iters_maxL.head(countL) );
  }
  if(warningL & warningK) {
    warnings = Rcpp::List::create(Rcpp::Named("top_maxK") = iters_maxK.head(countK),
                                  Rcpp::Named("top_maxL") = iters_maxL.head(countL) );
  }
  
  return Rcpp::List::create(Rcpp::Named("mu") = mu.t(),
                            Rcpp::Named("sigma2") = sigma2.t(),
                            Rcpp::Named("obs_cluster") = out_M.t(),
                            Rcpp::Named("distr_cluster") = out_S.t(),
                            Rcpp::Named("pi") = out_pi.t(),
                            Rcpp::Named("omega") = out_omega,
                            Rcpp::Named("alpha") = out_alpha,
                            Rcpp::Named("beta") = out_beta,
                            Rcpp::Named("maxK") = out_maxK,
                            Rcpp::Named("maxL") = out_maxL,
                            Rcpp::Named("warnings") = warnings);
}

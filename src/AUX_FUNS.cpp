// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "AUX_FUNS.h"
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

// sample from Dirichlet distribution
arma::vec rdirichlet(arma::vec par) 
{
  int distribution_size = par.n_elem ;
  arma::vec out = arma::zeros(distribution_size) ;
  double sum_term = 0 ;
  
  // draw Gamma variables
  for (int j = 0; j < distribution_size; j++) {
    double gam = R::rgamma(par[j], 1.0) ;
    out(j) = gam ;
    sum_term += gam ;
  }
  // normalize
  for (int j = 0; j < distribution_size; j++) { 
    out(j) = out(j)/sum_term ; 
  }
  return(out) ;
}

//---------------------------------------------------------

arma::vec relabel_arma(arma::vec cluster) 
{
  int n = cluster.n_elem ;
  arma::vec newlabels(n) ;
  int lab ;
  
  arma::vec uniquecl = arma::unique(cluster) ;
  int n_distinct = uniquecl.n_elem ;
  int maxcl = cluster.max() ;
  
  if( maxcl > (n_distinct - 1) ) {
    for(int k = 0; k < n_distinct ; k++) {
      lab = uniquecl[k] ;
      arma::uvec indlab = find(cluster == lab) ;
      for(unsigned int h = 0; h < indlab.n_elem; h++) { newlabels[indlab[h]] = k ; }
    }
  } else {
    newlabels = cluster ;
  }
  return(newlabels);
}


//---------------------------------------------------------

// sample fun
int sample_i(arma::vec ids, arma::vec prob)
{
  int out ;
  out = Rcpp::RcppArmadillo::sample(ids, 1, false, prob)[0] ;
  return(out) ;
}

//---------------------------------------------------------


// update model parameters <--- change if different likelihood
Rcpp::List sample_model_parameters(const arma::vec& y, 
                                  arma::vec M_iter,
                                  int maxL,
                                  double m0, double tau0,
                                  double lambda0, double gamma0)
{
  arma::vec out_mu(maxL) ;
  arma::vec out_sigma2(maxL) ;
  
  double new_m0 ; double new_tau0 ;
  double new_lambda0; double new_gamma0 ;
  
  for(int l = 0; l < maxL ; l++)
  {
    arma::uvec M_l = find( M_iter == l ) ;
    int n_l = M_l.n_elem ;
    
    if( n_l > 0) { 
      
      double sum_l = arma::accu( y(M_l) ); // sum of y_i

      // new parameters
      new_tau0 = tau0 + n_l ;
      new_m0 = (tau0 * m0 + sum_l)/new_tau0 ;
      new_lambda0 = lambda0 + n_l / 2 ;
      new_gamma0 = gamma0 + 0.5 * ( (n_l-1)*arma::var(y(M_l)) + (tau0*n_l)/new_tau0 * (sum_l/n_l - m0) * (sum_l/n_l - m0)  ) ;
      

      out_sigma2(l) = 1/R::rgamma( new_lambda0, 1/new_gamma0) ;
      out_mu(l) = R::rnorm( new_m0, sqrt(out_sigma2(l)/new_tau0) ) ;
      
    } else {
      out_sigma2(l) = 1/R::rgamma(lambda0, 1/gamma0) ;
      out_mu(l) = R::rnorm( m0, sqrt(out_sigma2(l)/tau0) )  ;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("out_mu") = out_mu,
                            Rcpp::Named("out_sigma2") = out_sigma2);
}

// Compute slice coefficients
// We use a deterministic sequence as in Denti et al. (2021)
// See Sec. C1 of the Supplementary Material for details
double fun_xi(double kappa, int i)
{
  double out = log(1-kappa) + (i-1) * log(kappa) ;
  return(exp(out)) ;
}

// Truncating threshold of the slice sampler
// See Sec. C1 of the Supplementary Material of Denti et al. (2021) for details
int compute_trunc(double u_min, double kappa)
{
  int out ;
  out = floor( (log(u_min) - log(1 - kappa)) / log(kappa) ) ;
  return(out) ;
}


// Compute stick-breaking weights starting from the vector of beta r.v. Beta(a_k,b_k)
// Sethuraman (1994) construction
arma::vec stick_breaking(arma::vec beta_var)
{
  int len = beta_var.n_elem ;
  arma::vec out(len) ;
  out(0) = beta_var(0) ;
  arma::vec cs_log_one_m_beta = arma::cumsum(log(1.0-beta_var));
  for(int k = 1; k < len; k++)
  {
    out(k) = exp( log(beta_var(k)) + cs_log_one_m_beta(k-1)) ;
  }
  return(out) ;
}

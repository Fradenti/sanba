#include <RcppArmadillo.h>
#include "D_CAVI_UPDATES.h"
#include "E_ELBO_COMPONENTS.h"

// [[Rcpp::export]]
Rcpp::List main_vb_fSAN_cpp(arma::field<arma::colvec> Y_grouped,
                             int const L,
                             int const K,
                             int const J,
                             arma::field<arma::mat> XI_ijl,
                             arma::mat RHO_jk,
                             arma::colvec Nj,
                             double m0,
                             double k0,
                             double a0,
                             double b0,
                             arma::colvec ml,
                             arma::colvec kl,
                             arma::colvec al,
                             arma::colvec bl,
                             double alpha_bar,
                             double beta_bar,
                             double epsilon,
                             int maxSIM,
                             bool verbose = 0){

  // random init
  arma::mat beta_lk(L,K); beta_lk.fill(beta_bar);
  arma::colvec alpha_k(K); alpha_k.fill(alpha_bar);

  arma::colvec  ELBO_val(maxSIM);

  arma::colvec alpha_star_k(K);
  arma::mat beta_star_lk(L,K);
  arma::mat var_par_theta(L,4);
  arma::mat var_par_v(K,3);

  int Q = maxSIM;
  //arma::mat XXX(maxSIM,6);

  for(int ii = 0; ii<maxSIM; ii++){
    R_CheckUserInterrupt();
    //Rcout << "*\n";
    beta_star_lk = Update_beta_dirlk_cpp(XI_ijl,
                                         RHO_jk,
                                         beta_lk,
                                         L,
                                         J,
                                         K);


    // M
    XI_ijl = Update_XIijl_cpp_fiSAN(Y_grouped,
                                    RHO_jk,
                                    beta_star_lk,
                                    Nj,
                                    ml,
                                    kl,
                                    al,
                                    bl,
                                    L,
                                    J,
                                    K);
    //Rcout << "3.5*\n";

    //Rcout << "*\n";

    // THETA
    var_par_theta  = Update_THETAl_cpp(Y_grouped,
                                       XI_ijl,
                                       m0,
                                       k0,
                                       a0,
                                       b0,
                                       L,
                                       J);
    ml = var_par_theta.col(0);
    kl = var_par_theta.col(1);
    al = var_par_theta.col(2);
    bl = var_par_theta.col(3);


    // V
    alpha_star_k = Update_alpha_dirk_cpp(RHO_jk, // collection of  J objectes: nj*L matrices
                                         alpha_k);


    // S
    RHO_jk = Update_RHOjk_cpp_overCAM(XI_ijl, // collection of  J objectes: nj*L matrices
                                    alpha_star_k,
                                    beta_star_lk,
                                    L,
                                    J,
                                    K);

    double ELBO_OMEGA = elbo_p_omega(beta_star_lk,
                                     beta_lk, L, K) -
                                       elbo_q_omega(beta_star_lk,
                                                    L,
                                                    K);

    double ELBO_M = elbo_p_M_fiSAN(XI_ijl,
                                   RHO_jk,
                                   beta_star_lk,
                                   L, K, J) -
                                     elbo_q_M(XI_ijl, J);

    double ELBO_PI = elbo_p_pi(alpha_star_k,  alpha_k) -
      elbo_q_pi(alpha_star_k);

    double ELBO_S = elbo_p_S_overCAM(RHO_jk,alpha_star_k) - elbo_q_S(RHO_jk);

    double  Elbo_pLIK = elbo_p_Y(Y_grouped,
                                 XI_ijl,
                                 ml,kl,al,bl,
                                 L,J);
    double ELBO_THETA  = elbo_p_THETA(m0, k0, a0, b0,
                                      ml, kl, al, bl)-
                                        elbo_q_THETA(ml, kl,al,bl);

    ELBO_val[ii] =  Elbo_pLIK +
      ELBO_S + ELBO_M + ELBO_PI + ELBO_OMEGA + ELBO_THETA;


    if(ii>1) {
      double diff = (ELBO_val[ii]-ELBO_val[ii-1]);
      if(verbose){
        Rcpp::Rcout << "Iteration #" << ii << " - Elbo increment: " << diff << "\n";
      }
    if(((diff<0) & (std::fabs(diff) > 1e-5))){
        Rcpp::Rcout << "Warning! Iteration #" << ii << " presents an ELBO decrement!\n";
      }
      if( diff < epsilon ) {
        if(verbose){
          Rcpp::Rcout << "Final Elbo value: " << (ELBO_val[ii]) << "\n";
        }
        Q = ii;
        break;
      }
    }

  }

  if(Q == maxSIM){Q = Q-1;}

  arma::colvec ELBO_v2 = ELBO_val.rows(1,Q);
  //arma::mat YYY = XXX.rows(1,Q);

  Rcpp::List results = Rcpp::List::create(
    Rcpp::_["theta_l"]  = var_par_theta,
    Rcpp::_["XI"] = XI_ijl,
    Rcpp::_["RHO"] = RHO_jk,
    Rcpp::_["a_dirichlet_k"] = alpha_star_k,
    Rcpp::_["b_dirichlet_lk"] = beta_star_lk,
    Rcpp::_["Elbo_val"] = ELBO_v2
  );
  return(results);
}


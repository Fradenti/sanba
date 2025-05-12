#ifndef E_ELBO_COMPONENTS
#define E_ELBO_COMPONENTS

#include <RcppArmadillo.h>
#include "B_EXP_VALUES.h"
#include "C_NORM_CONST.h"


double elbo_p_M_fiSAN(arma::field<arma::mat> XI_ijl,
                      arma::mat RHO_jk,
                      arma::mat beta_star_lk,
                      int const L,
                      int const K,
                      int const J);


double elbo_p_M_CAM(arma::field<arma::mat> XI_ijl,
                    arma::mat RHO_jk,
                    arma::mat ElnOM_lk,
                    int const L,
                    int const K,
                    int const J);


double elbo_p_S_overCAM(arma::mat RHO_jk,
                        arma::colvec alpha_star_k);
// Common

double elbo_p_Y(arma::field<arma::colvec> Y_grouped,
                arma::field<arma::mat> XI_ijl,
                arma::colvec ml,
                arma::colvec kl,
                arma::colvec al,
                arma::colvec bl,
                int L,
                int J);

// -----------------------------------------------------------------------------

double elbo_q_M(arma::field<arma::mat> XI_ijl,
                int J);

// -----------------------------------------------------------------------------

double elbo_p_v(arma::colvec a_tilde_k,
                arma::colvec b_tilde_k,
                double const a_tilde,
                double const b_tilde,
                int const K);

// -----------------------------------------------------------------------------

double elbo_q_v(arma::colvec a_tilde_k,
                arma::colvec b_tilde_k,
                int const K);

// -----------------------------------------------------------------------------

double elbo_p_THETA(double m0,
                    double k0,
                    double a0,
                    double b0,
                    arma::colvec ml,
                    arma::colvec kl,
                    arma::colvec al,
                    arma::colvec bl);


// -----------------------------------------------------------------------------

double elbo_q_THETA(arma::colvec ml,
                    arma::colvec kl,
                    arma::colvec al,
                    arma::colvec bl);

// -----------------------------------------------------------------------------

double elbo_p_S(arma::mat RHO_jk,
                arma::colvec ElnPI);

// -----------------------------------------------------------------------------

double elbo_q_S(arma::mat RHO_jk);

// -----------------------------------------------------------------------------

double elbo_p_omega(arma::mat beta_star_lk,
                    arma::mat beta_lk,
                    int const L,
                    int const K);

// -----------------------------------------------------------------------------

double elbo_q_omega(arma::mat beta_star_lk,
                    int const L,
                    int const K);

// -----------------------------------------------------------------------------

double elbo_p_U(arma::mat a_bar_Ulk,
                arma::mat b_bar_Ulk,
                double const a_bar,
                double const b_bar,
                int const L,
                int const K);

// -----------------------------------------------------------------------------

double elbo_q_U(arma::mat a_bar_Ulk,
                arma::mat b_bar_Ulk,
                int const L,
                int const K);

// -----------------------------------------------------------------------------

double elbo_p_pi(arma::colvec alpha_star_k,
                 arma::colvec alpha_k);

// -----------------------------------------------------------------------------

double elbo_q_pi(arma::colvec alpha_star_k);

// -----------------------------------------------------------------------------

double elbo_conc_par_CAM(arma::colvec conc_hyper,
                     arma::colvec S_concDP);

// -----------------------------------------------------------------------------

double elbo_conc_par_fiSAN(arma::colvec conc_hyper,
                           arma::colvec S_concDP);

// -----------------------------------------------------------------------------

double elbo_p_v_CP(arma::colvec a_tilde_k,
                   arma::colvec b_tilde_k,
                  // double const a_tilde,
                   double const b_tilde,
                   arma::colvec S_concDP,
                   int const K);

// -----------------------------------------------------------------------------

double elbo_p_U_CP(arma::mat a_bar_Ulk,
                   arma::mat b_bar_Ulk,
                   //double const a_bar,
                   double const b_bar,
                   arma::colvec S_concDP,
                   int const L, int const K);


#endif

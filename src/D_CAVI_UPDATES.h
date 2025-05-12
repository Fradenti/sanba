#ifndef D_CAVI_UPDATES
#define D_CAVI_UPDATES

#include <RcppArmadillo.h>

#include "A_AUX.h"

#include "B_EXP_VALUES.h"

arma::mat Update_Vk_cpp(int const K,
                        double a_tilde,
                        double b_tilde,
                        arma::mat RHO_jk);

// -----------------------------------------------------------------------------

arma::mat Update_THETAl_cpp(arma::field<arma::colvec> Y_grouped,
                            arma::field<arma::mat> XI_ijl,
                            double m0,
                            double k0,
                            double a0,
                            double b0,
                            int const L,
                            int const J);

// -----------------------------------------------------------------------------

arma::mat Update_beta_dirlk_cpp(arma::field<arma::mat> XI_ijl,
                                arma::mat RHO_jk,
                                arma::mat beta_lk,
                                int const L,
                                int const J,
                                int const K);

// -----------------------------------------------------------------------------

arma::colvec Update_alpha_dirk_cpp(arma::mat RHO_jk,
                                arma::colvec alpha_k);

// -----------------------------------------------------------------------------

arma::cube Update_Ulk_cpp(arma::field<arma::mat> XI_ijl,
                          arma::mat RHO_jk,
                          double const a_bar,
                          double const b_bar,
                          int const L,
                          int const J,
                          int const K);



/// fiSAN - specific

arma::mat Update_RHOjk_cpp_fiSAN(arma::field<arma::mat> XI_ijl, // collection of  J objectes: nj*L matrices
                                 arma::colvec ElnPI_k,
                                 arma::mat beta_star_lk,
                                 int const L,
                                 int const J,
                                 int const K);

// -----------------------------------------------------------------------------

arma::field<arma::mat> Update_XIijl_cpp_fiSAN(arma::field<arma::colvec> Y_grouped,
                                              arma::mat RHO_jk,
                                              arma::mat beta_star_lk,
                                              arma::colvec Nj,
                                              arma::colvec ml,
                                              arma::colvec kl,
                                              arma::colvec al,
                                              arma::colvec bl,
                                              int const L,
                                              int const J,
                                              int const K);


/// CAM - specific
arma::mat  Update_RHOjk_cpp_CAM(arma::field<arma::mat> XI_ijl,
                                arma::colvec ElnPI_k,
                                arma::mat    ElnOM_lk,
                                int const L,
                                int const J,
                                int const K);

// -----------------------------------------------------------------------------

arma::field<arma::mat> Update_XIijl_cpp_CAM(arma::field<arma::colvec> Y_grouped,
                                            arma::mat RHO_jk,
                                            arma::mat ElnOM_lk,
                                            arma::colvec Nj,
                                            arma::colvec ml,
                                            arma::colvec kl,
                                            arma::colvec al,
                                            arma::colvec bl,
                                            int const L,
                                            int const J,
                                            int const K);

// -----------------------------------------------------------------------------

arma::mat Update_RHOjk_cpp_overCAM(arma::field<arma::mat> XI_ijl, // collection of  J objectes: nj*L matrices
                                   arma::colvec alpha_star_k,
                                   arma::mat beta_star_lk,
                                   int const L,
                                   int const J,
                                   int const K);

// -----------------------------------------------------------------------------

arma::colvec Update_s_concentration_par(arma::colvec a_tilde_Vk,
                                        arma::colvec b_tilde_Vk,
                                        arma::mat a_bar_Ulk,
                                        arma::mat b_bar_Ulk,
                                        arma::colvec conc_hyper,
                                        int L,
                                        int K);

// -----------------------------------------------------------------------------

arma::colvec Update_s_concentration_par_fiSAN(arma::colvec a_tilde_Vk,
                                              arma::colvec b_tilde_Vk,
                                              arma::colvec conc_hyper,
                                              int K);
#endif


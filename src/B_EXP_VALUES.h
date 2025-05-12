#ifndef B_EXP_VALUES
#define B_EXP_VALUES

#include <RcppArmadillo.h>

arma::colvec E_log_beta(arma:: colvec a,
                        arma:: colvec b);

arma::colvec E_log_IG(arma::colvec a,
                      arma::colvec b);
  
arma::colvec E_log_DIR(arma::colvec a);
  

arma::colvec E_log_p_Y_Mtheta_cpp(arma::colvec Y,
                                  double ml,
                                  double kl,
                                  double al,
                                  double bl);

#endif


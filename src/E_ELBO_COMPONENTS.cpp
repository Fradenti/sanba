#include "E_ELBO_COMPONENTS.h"

// fiSAN specific - used by overcam as well

double elbo_p_M_fiSAN(arma::field<arma::mat> XI_ijl,
                      arma::mat RHO_jk,
                      arma::mat beta_star_lk,
                      int const L,
                      int const K,
                      int const J){

  arma::mat N_jl(J,L);
  arma::mat E_ln_omega(L,K);

  for(int k=0; k<K; k++){
    E_ln_omega.col(k) =
      E_log_DIR(beta_star_lk.col(k)); // L x K
  }
  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_ijl(j),0); //J x L
  }

  arma::mat X_jl  = RHO_jk * E_ln_omega.t();
  double Z = accu( N_jl % X_jl);
  return(Z);
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// CAM specific

double elbo_p_M_CAM(arma::field<arma::mat> XI_ijl,
                    arma::mat RHO_jk,
                    arma::mat ElnOM_lk,
                    int const L,
                    int const K,
                    int const J){

  arma::mat N_jl(J,L);

  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_ijl(j),0); //J x L
  }

  arma::mat X_jl  = RHO_jk * ElnOM_lk.t();
  double Z = accu( N_jl % X_jl);
  return(Z);
}

// overcam

double elbo_p_S_overCAM(arma::mat RHO_jk,
                        arma::colvec alpha_star_k){


  arma::colvec ElnPi = E_log_DIR(alpha_star_k);
  double Z  = arma::accu( RHO_jk * ElnPi );
  return(Z);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Common

double elbo_p_Y(arma::field<arma::colvec> Y_grouped,
                arma::field<arma::mat> XI_ijl,
                arma::colvec ml,
                arma::colvec kl,
                arma::colvec al,
                arma::colvec bl,
                int L,
                int J){

  double Z = 0;
  arma::colvec Q = E_log_IG(al,bl);

  for(int j=0; j < J; j++){
    for(int l=0; l < L; l++){

      arma::colvec v =   (Q[l] + 1.0/kl[l] + //Nj x 1
        ( Y_grouped[j] - ml[l]) % ( Y_grouped[j] - ml[l] ) * ( al[l]/bl[l] ) );
      Z += arma::accu(XI_ijl[j].col(l) % v);

    }
  }
  return(-.5*Z);

}

// -----------------------------------------------------------------------------

double elbo_q_M(arma::field<arma::mat> XI_ijl, int J){


  arma::colvec Z(J);
  for(int j=0; j<J; j++){
    Z(j) = arma::accu(XI_ijl[j] % log(XI_ijl[j] + 1e-12) );
  }
  double X = arma::accu(Z);

  return(X);
}

// -----------------------------------------------------------------------------

double elbo_p_v(arma::colvec a_tilde_k,
                arma::colvec b_tilde_k,
                double const a_tilde,
                double const b_tilde, int const K){

  a_tilde_k.shed_row(K-1);
  b_tilde_k.shed_row(K-1);

  arma::colvec Y = lbeta_normconst_cpp(a_tilde,b_tilde) +
    E_log_beta(a_tilde_k,b_tilde_k) * ( a_tilde - 1 )+
    E_log_beta(b_tilde_k,a_tilde_k) * ( b_tilde - 1 );

  return(arma::accu(Y));
}


double elbo_p_v_CP(arma::colvec a_tilde_k,
                arma::colvec b_tilde_k,
                //double const a_tilde,
                double const b_tilde,
                arma::colvec S_concDP,
                int const K){

  a_tilde_k.shed_row(K-1);
  b_tilde_k.shed_row(K-1);

  arma::colvec Y =
    // E_log_beta(a_tilde_k,b_tilde_k) * ( a_tilde - 1 )+
    E_log_beta(b_tilde_k,a_tilde_k) * ( b_tilde - 1 );

  return((K-1)  * ( R::digamma(S_concDP[0]) - log(S_concDP[1]) ) + arma::accu(Y));
}



// -----------------------------------------------------------------------------

double elbo_q_v(arma::colvec a_tilde_k,
                arma::colvec b_tilde_k, int const K){
  a_tilde_k.shed_row(K-1);
  b_tilde_k.shed_row(K-1);


  arma::colvec Y = lbeta_normconst_vec_cpp(a_tilde_k,b_tilde_k) +
    E_log_beta(a_tilde_k,b_tilde_k) % (a_tilde_k - 1 )+
    E_log_beta(b_tilde_k,a_tilde_k) % (b_tilde_k - 1 );

  return(arma::accu(Y));
}

// -----------------------------------------------------------------------------

double elbo_p_THETA(double m0, double k0, double a0, double b0,
                    arma::colvec ml, arma::colvec kl,
                    arma::colvec al, arma::colvec bl){

  double  Z = - (a0 + 1.5) * arma::accu( E_log_IG(al,bl) ) -
    b0        * arma::accu( al / bl ) -
    k0 / 2    * arma::accu(1 / kl + (al / bl) % (ml-m0) % (ml-m0));

  return(Z);
}


// -----------------------------------------------------------------------------

double elbo_q_THETA(arma::colvec ml, arma::colvec kl,
                    arma::colvec al, arma::colvec bl){

  double K = arma::accu(al % log(bl) - lgamma(al) +.5* log(kl) -
                        (al + 1.5) % E_log_IG(al,bl) - al );

  return(K);
}

// -----------------------------------------------------------------------------

double elbo_p_S(arma::mat RHO_jk,
                arma::colvec ElnPI){

  arma::colvec mdot_k = arma::sum(RHO_jk, 0).t();
  double Z = arma::accu(mdot_k % ElnPI);
  return(Z);
}

// -----------------------------------------------------------------------------

double elbo_q_S(arma::mat RHO_jk){
  double y = arma::accu(RHO_jk % log(RHO_jk + 1e-12));
  return(y);
}

// -----------------------------------------------------------------------------

double elbo_p_omega(arma::mat beta_star_lk,
                    arma::mat beta_lk,
                    int const L, int const K){

  arma::mat E_ln_omega(L,K);
  double Konsts = arma::accu(ldirichlet_normconst_vec_cpp(beta_lk));
  for(int k=0; k<K; k++){
    E_ln_omega.col(k) =
      E_log_DIR(beta_star_lk.col(k)); // L x K
  }

  double G =  Konsts + arma::accu( E_ln_omega % (beta_lk - 1.0) );
  return(G);
}

// -----------------------------------------------------------------------------

double elbo_q_omega(arma::mat beta_star_lk,
                    int const L, int const K){

  arma::mat E_ln_omega(L,K);
  double Konsts = arma::accu(ldirichlet_normconst_vec_cpp(beta_star_lk));

  for(int k=0; k<K; k++){
    E_ln_omega.col(k) =
      E_log_DIR(beta_star_lk.col(k)); // L x K
  }

  double G =  Konsts + arma::accu( E_ln_omega % (beta_star_lk - 1.0) );
  return(G);
}

// -----------------------------------------------------------------------------

double elbo_p_U(arma::mat a_bar_Ulk,
                arma::mat b_bar_Ulk,
                double const a_bar, double const b_bar,
                int const L, int const K){

  a_bar_Ulk.shed_row(L-1);
  b_bar_Ulk.shed_row(L-1);

  arma::colvec R(K);

  for(int k = 0; k < K; k ++){

    R[k] =  (a_bar-1) * arma::accu( E_log_beta(a_bar_Ulk.col(k),b_bar_Ulk.col(k))) +
            (b_bar-1) * arma::accu( E_log_beta(b_bar_Ulk.col(k),a_bar_Ulk.col(k)));

  }

  return(arma::accu(R));

}


double elbo_p_U_CP(arma::mat a_bar_Ulk,
                   arma::mat b_bar_Ulk,
                   //double const a_bar,
                   double const b_bar,
                   arma::colvec S_concDP,
                   int const L, int const K){

  a_bar_Ulk.shed_row(L-1);
  b_bar_Ulk.shed_row(L-1);

  arma::colvec R(K);
  double C = (L-1) * K * ( R::digamma(S_concDP[2]) - log(S_concDP[3]) );

  for(int k = 0; k < K; k ++){

    R[k] =  //arma::accu( (a_bar-1) * E_log_beta(a_bar_Ulk.col(k),b_bar_Ulk.col(k))) +
              arma::accu(  E_log_beta(b_bar_Ulk.col(k),a_bar_Ulk.col(k)));

  }

  return(C + (b_bar - 1) * arma::accu(R));

}




// -----------------------------------------------------------------------------

double elbo_q_U(arma::mat a_bar_Ulk,
                arma::mat b_bar_Ulk,
                int const L, int const K){

  a_bar_Ulk.shed_row(L-1);
  b_bar_Ulk.shed_row(L-1);
  arma::colvec R(K);
  double C  = arma::accu(lbeta_normconst_mat_cpp(a_bar_Ulk,b_bar_Ulk));

  for(int k = 0; k < K; k ++){

    R[k] =  arma::accu( (a_bar_Ulk.col(k)-1) % E_log_beta(a_bar_Ulk.col(k),b_bar_Ulk.col(k))) +
            arma::accu( (b_bar_Ulk.col(k)-1) % E_log_beta(b_bar_Ulk.col(k),a_bar_Ulk.col(k)));

  }

  return(C + arma::accu(R));

}

// -----------------------------------------------------------------------------

double elbo_p_pi(arma::colvec alpha_star_k,
                 arma::colvec alpha_k){

  double Konsts = ldirichlet_normconst_cpp(alpha_k);
  double G =  Konsts + arma::accu( E_log_DIR(alpha_star_k) % (alpha_k - 1) );
  return(G);
}

// -----------------------------------------------------------------------------

double elbo_q_pi(arma::colvec alpha_star_k){

  double Konsts = ldirichlet_normconst_cpp(alpha_star_k);

  double G =  Konsts + arma::accu( E_log_DIR(alpha_star_k) % (alpha_star_k - 1) );
  return(G);
}

// -----------------------------------------------------------------------------

double elbo_conc_par_CAM(arma::colvec conc_hyper,
                     arma::colvec S_concDP){
  double pa =
  // p(gamma_alpha)
  conc_hyper[0] * log(conc_hyper[1]) - lgamma(conc_hyper[0]) +
  (conc_hyper[0]-1) * (R::digamma(S_concDP[0]) - log(S_concDP[1])) -
   conc_hyper[1] * S_concDP[0]/S_concDP[1];
  // q(gamma_alpha)
  double qa = S_concDP[0] * log(S_concDP[1]) - lgamma(S_concDP[0]) +
  (S_concDP[0]-1) * (R::digamma(S_concDP[0])-log(S_concDP[1])) -  S_concDP[0];


  // p(gamma_beta)
  double pb =
    conc_hyper[2] * log(conc_hyper[3]) - lgamma(conc_hyper[2]) +
  (conc_hyper[2]-1) * (R::digamma(S_concDP[2])-log(S_concDP[3])) -
  conc_hyper[3] * S_concDP[2]/S_concDP[3];
  // q(gamma_beta)
  double qb = S_concDP[2] * log(S_concDP[3]) - lgamma(S_concDP[2]) +
  (S_concDP[2]-1) * (R::digamma(S_concDP[2])-log(S_concDP[3])) -  S_concDP[2] ;

  //Rcpp::Rcout << pa << " -- " << pb << " -- " << qa << " -- " << qb << " -- \n" ;
  return(pa + pb - qa - qb);
}

// -----------------------------------------------------------------------------

double elbo_conc_par_fiSAN(arma::colvec conc_hyper,
                         arma::colvec S_concDP){
  double pa =
    // p(gamma_alpha)
    conc_hyper[0] * log(conc_hyper[1]) - lgamma(conc_hyper[0]) +
    (conc_hyper[0]-1) * (R::digamma(S_concDP[0]) - log(S_concDP[1])) -
    conc_hyper[1] * S_concDP[0]/S_concDP[1];
  // q(gamma_alpha)
  double qa = S_concDP[0] * log(S_concDP[1]) - lgamma(S_concDP[0]) +
    (S_concDP[0]-1) * (R::digamma(S_concDP[0])-log(S_concDP[1])) -  S_concDP[0];

  return(pa - qa);
}

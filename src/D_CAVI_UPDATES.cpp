#include "D_CAVI_UPDATES.h"

arma::mat Update_Vk_cpp(int const K,
                        double  a_tilde,
                        double  b_tilde,
                        arma::mat RHO_jk){

  arma::colvec mk = (arma::sum(RHO_jk,0)).t();
  arma::colvec rev_cs_mk = reverse_cumsum_cpp(mk);

  arma::colvec a_tilde_vk      = mk        + a_tilde;
  a_tilde_vk[K-1] = 1;
  arma::colvec b_tilde_vk      = rev_cs_mk + b_tilde;
  b_tilde_vk[K-1] = 1e-10;

  arma::colvec E_ln_Vk    = E_log_beta(a_tilde_vk, b_tilde_vk);
  arma::colvec E_ln_1mVk  = E_log_beta(b_tilde_vk, a_tilde_vk);
  arma::colvec sE_ln_1mVk = shift(E_ln_1mVk, +1);

  sE_ln_1mVk[0] = 0;

  arma::colvec CS_E_ln_1mVk = arma::cumsum(sE_ln_1mVk);
  arma::mat results(K,3);

  results.col(0) = a_tilde_vk;
  results.col(1) = b_tilde_vk;
  results.col(2) = E_ln_Vk + CS_E_ln_1mVk ;

  return(results);
}

/* maybe not needed for now
arma::mat Update_Vk_DP_conc_cpp(int const K,
                                double const s1_star,
                                double const s2_star,
                                arma::mat RHO_jk){

  arma::colvec mk = (arma::sum(RHO_jk,0)).t();
  arma::colvec rev_cs_mk = reverse_cumsum_cpp(mk);

  arma::colvec a_tilde_vk      = mk        + 1;
  a_tilde_vk[K-1] = 1;
  arma::colvec b_tilde_vk      = rev_cs_mk + s1_star/s2_star;
  b_tilde_vk[K-1] = 1e-10;

  arma::colvec E_ln_Vk    = E_log_beta(a_tilde_vk, b_tilde_vk);
  arma::colvec E_ln_1mVk  = E_log_beta(b_tilde_vk, a_tilde_vk);
  arma::colvec sE_ln_1mVk = shift(E_ln_1mVk, +1);

  sE_ln_1mVk[0] = 0;

  arma::colvec CS_E_ln_1mVk = arma::cumsum(sE_ln_1mVk);
  arma::mat results(K,3);

  results.col(0) = a_tilde_vk;
  results.col(1) = b_tilde_vk;
  results.col(2) = E_ln_Vk + CS_E_ln_1mVk ;

  return(results);
}
*/


// -----------------------------------------------------------------------------

arma::mat Update_THETAl_cpp(arma::field<arma::colvec> Y_grouped,
                            arma::field<arma::mat> XI_ijl,
                            double m0, double k0,
                            double a0, double b0,
                            int const L,int const J){

  arma::mat results(L,4);

  arma::rowvec Nl(L);
  arma::rowvec SumY_l(L);
  arma::rowvec SumY2_l(L);

  Nl.zeros();
  SumY_l.zeros();
  SumY2_l.zeros();

  for(int j=0; j <J; j++){
    arma::colvec subY = Y_grouped[j];
    Nl      += arma::sum( XI_ijl(j), 0);
    SumY_l  += subY.t() * XI_ijl(j);
    SumY2_l += pow(subY.t(),2) * XI_ijl(j);
  }

  for(int l = 0; l<L; l++){
      double Ybar_l = 0;
      double S2_l   = 0;

        if(Nl[l] > 0){
          Ybar_l = SumY_l[l]/Nl[l];
          S2_l   = SumY2_l[l] - pow(Ybar_l,2) * Nl[l];
        }

          double    kl = (k0 + Nl[l]);
          results(l,0) = (m0 * k0 + SumY_l[l]) / (kl) ;
          results(l,1) = kl ;
          results(l,2) = a0 + Nl[l] / 2 ;
          results(l,3) = b0 + 0.5 * (S2_l +
                                       ( ( (k0 * Nl[l]) / kl ) * pow(Ybar_l - m0,2)) ) ;

    }

  return(results);

}

// -----------------------------------------------------------------------------

arma::mat Update_beta_dirlk_cpp(arma::field<arma::mat> XI_ijl,
                                arma::mat RHO_jk,
                                arma::mat beta_lk,
                                int const L,
                                int const J,
                                int const K){

  arma::mat N_lj(L,J);
  arma::mat beta_star_lk(L,K);

  for(int j=0; j<J; j++){
    N_lj.col(j) = arma::sum(XI_ijl(j),0).t();
  }

  arma::mat Q_lk = N_lj * RHO_jk ;

  for(int k = 0; k<K; k++){

    beta_star_lk.col(k) = beta_lk.col(k) + Q_lk.col(k);

  }

  return(beta_star_lk);
}



arma::colvec Update_alpha_dirk_cpp(arma::mat RHO_jk,
                                   arma::colvec alpha_k){

  arma::rowvec Nk = arma::sum(RHO_jk,0);

  return(alpha_k + Nk.t());
}

// -----------------------------------------------------------------------------

arma::cube Update_Ulk_cpp(arma::field<arma::mat> XI_ijl,
                          arma::mat RHO_jk,
                          double const a_bar,
                          double const b_bar,
                          int const L,
                          int const J,
                          int const K){
  //to improve
  arma::mat N_jl(J,L);


  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_ijl(j),0);
  }

  arma::mat Q_lk = N_jl.t() * RHO_jk ;
  arma::mat a_bar_Ulk(L,K);
  arma::mat b_bar_Ulk(L,K);
  arma::mat E_lnOmega_lk(L,K);


  for(int k = 0; k<K; k++){

    arma::colvec rQl_k  = reverse_cumsum_cpp(Q_lk.col(k));
    arma::colvec G1 = a_bar + Q_lk.col(k);
    arma::colvec G2 = b_bar + rQl_k;

    G1[L-1] = 1;
    G2[L-1] = 1e-10;

    a_bar_Ulk.col(k) =  G1;
    b_bar_Ulk.col(k) =  G2;

    arma::colvec E_ln_Ul_k    = E_log_beta(G1, G2);
    arma::colvec E_ln_1mUl_k  = E_log_beta(G2, G1);
    arma::colvec sE_ln_1mUl_k = shift(E_ln_1mUl_k, +1);
    sE_ln_1mUl_k[0] = 0;

    arma::colvec CS_E_ln_1mUL_k = arma::cumsum(sE_ln_1mUl_k);
    E_lnOmega_lk.col(k) = CS_E_ln_1mUL_k + E_ln_Ul_k;
  }


  arma::cube results(L,K,3);
  results.slice(0) = a_bar_Ulk;
  results.slice(1) = b_bar_Ulk;
  results.slice(2) = E_lnOmega_lk;

  return(results);
}


/* maybe not needed for now
arma::cube Update_Ulk_DP_conc_cpp(arma::field<arma::mat> XI_ijl,
                                   arma::mat RHO_jk,
                                   double const s3_star,
                                   double const s4_star,
                                   int const L,
                                   int const J,
                                   int const K){
  //to improve
  arma::mat N_jl(J,L);


  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_ijl(j),0);
  }

  arma::mat Q_lk = N_jl.t() * RHO_jk ;
  arma::mat a_bar_Ulk(L,K);
  arma::mat b_bar_Ulk(L,K);
  arma::mat E_lnOmega_lk(L,K);


  for(int k = 0; k<K; k++){

    arma::colvec rQl_k  = reverse_cumsum_cpp(Q_lk.col(k));
    arma::colvec G1 = 1 + Q_lk.col(k);
    arma::colvec G2 = s3_star/s4_star + rQl_k;

    G1[L-1] = 1;
    G2[L-1] = 1e-10;

    a_bar_Ulk.col(k) =  G1;
    b_bar_Ulk.col(k) =  G2;

    arma::colvec E_ln_Ul_k    = E_log_beta(G1, G2);
    arma::colvec E_ln_1mUl_k  = E_log_beta(G2, G1);
    arma::colvec sE_ln_1mUl_k = shift(E_ln_1mUl_k, +1);
    sE_ln_1mUl_k[0] = 0;

    arma::colvec CS_E_ln_1mUL_k = arma::cumsum(sE_ln_1mUl_k);
    E_lnOmega_lk.col(k) = CS_E_ln_1mUL_k + E_ln_Ul_k;
  }


  arma::cube results(L,K,3);
  results.slice(0) = a_bar_Ulk;
  results.slice(1) = b_bar_Ulk;
  results.slice(2) = E_lnOmega_lk;

  return(results);
}
*/




// fiSAN - specific
arma::mat Update_RHOjk_cpp_fiSAN(arma::field<arma::mat> XI_ijl, // collection of  J objectes: nj*L matrices
                                 arma::colvec ElnPI_k,
                                 arma::mat beta_star_lk,
                                 int const L,
                                 int const J,
                                 int const K){

  arma::mat N_jl(J,L);


  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_ijl(j),0);
  }

  arma::mat unn_log_RHO_jk(J,K);
  arma::mat log_RHO_jk(J,K);
  arma::mat RHO_jk(J,K);
  arma::mat E_ln_omega(L,K);

  for(int k=0; k<K; k++){
    E_ln_omega.col(k) =
      E_log_DIR(beta_star_lk.col(k)); // L x K
  }

  arma::mat Z = N_jl * E_ln_omega;

  for(int k = 0; k < K; k++){
    unn_log_RHO_jk.col(k) =  Z.col(k) +  ElnPI_k(k);
  }
  for(int j = 0; j < J; j++){
    log_RHO_jk.row(j) = unn_log_RHO_jk.row(j) - LogSumExp_cpp(unn_log_RHO_jk.row(j));
  }

  return(exp(log_RHO_jk));
}


arma::colvec Update_s_concentration_par(arma::colvec a_tilde_Vk,
                                        arma::colvec b_tilde_Vk,
                                        arma::mat a_bar_Ulk,
                                        arma::mat b_bar_Ulk,
                                        arma::colvec conc_hyper,
                                        int L,
                                        int K){

  arma::colvec upd_par(4);
  a_bar_Ulk.shed_row(L-1);
  b_bar_Ulk.shed_row(L-1);
  a_tilde_Vk.shed_row(K-1);
  b_tilde_Vk.shed_row(K-1);

  arma::colvec R(K);

  for(int k = 0; k < K; k ++){

    R[k] =  arma::accu( E_log_beta(b_bar_Ulk.col(k),a_bar_Ulk.col(k)) );

  }

  upd_par[0] = conc_hyper[0] + K - 1 ; 
  upd_par[1] = conc_hyper[1] - arma::accu(E_log_beta(b_tilde_Vk,a_tilde_Vk));
  upd_par[2] = conc_hyper[2] + K * (L - 1);
  upd_par[3] = conc_hyper[3] - arma::accu(R);


  return(upd_par);
}

arma::colvec Update_s_concentration_par_fiSAN(arma::colvec a_tilde_Vk,
                                        arma::colvec b_tilde_Vk,
                                        arma::colvec conc_hyper,
                                        int K){

  arma::colvec upd_par(2);
  a_tilde_Vk.shed_row(K-1);
  b_tilde_Vk.shed_row(K-1);

  upd_par[0] = conc_hyper[0] + K - 1 ;
  upd_par[1] = conc_hyper[1] - arma::accu(E_log_beta(b_tilde_Vk,a_tilde_Vk));

  return(upd_par);
}


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
                                              int const K){


  arma::field<arma::mat> XI_ijl(J); // J different nj x L matrices
  arma::mat E_ln_omega(L,K);

  for(int k=0; k<K; k++){
    E_ln_omega.col(k) =
      E_log_DIR(beta_star_lk.col(k)); // L x K
  }

  arma::mat M_JL  =  RHO_jk * E_ln_omega.t();
  //to improve
  arma::rowvec logunn(L);

  for(int j= 0; j<J; j++){

    arma::mat TT(Nj[j],L);

    for(int l= 0; l<L; l++){
      TT.col(l) = E_log_p_Y_Mtheta_cpp(Y_grouped[j],
             ml[l], kl[l], al[l], bl[l]);
    }
    arma::mat temp1 = TT;  // Nj x L
    arma::mat tempres(Nj[j],L);

    for(int i=0; i<Nj[j]; i++){

      logunn = M_JL.row(j) + temp1.row(i);

      tempres.row(i) =  logunn - LogSumExp_cpp(logunn);
    }

    XI_ijl(j) = exp(tempres);

  }

  return(XI_ijl);
}

// CAM - specific
arma::mat  Update_RHOjk_cpp_CAM(arma::field<arma::mat> XI_ijl,
                               arma::colvec ElnPI_k,
                               arma::mat    ElnOM_lk,
                               int const L,
                               int const J,
                               int const K){


  arma::mat N_jl(J,L);


  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_ijl(j),0);
  }


  arma::mat unn_log_RHO_jk(J,K);
  arma::mat log_RHO_jk(J,K);
  arma::mat RHO_jk(J,K);

  arma::mat Z = N_jl * ElnOM_lk;

  for(int k = 0; k < K; k++){
    unn_log_RHO_jk.col(k) =  Z.col(k) +  ElnPI_k(k);
  }
  for(int j = 0; j < J; j++){
    log_RHO_jk.row(j) = unn_log_RHO_jk.row(j) - LogSumExp_cpp(unn_log_RHO_jk.row(j));
  }


  return(exp(log_RHO_jk));
}

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
                                        int const K){


  arma::field<arma::mat> XI_ijl(J); // J different nj x L matrices
  arma::mat M_LJ  =  ElnOM_lk * RHO_jk.t();
  arma::mat M_JL  = M_LJ.t();
  //to improve
  arma::rowvec logunn(L);

  for(int j= 0; j<J; j++){

    arma::mat TT(Nj[j],L);

    for(int l= 0; l<L; l++){
      TT.col(l) = E_log_p_Y_Mtheta_cpp(Y_grouped[j],
             ml[l], kl[l], al[l], bl[l]);
    }
    arma::mat temp1 = TT;  // Nj x L
    arma::mat tempres(Nj[j],L);

    for(int i=0; i<Nj[j]; i++){

      logunn = M_JL.row(j) + temp1.row(i);

      tempres.row(i) =  logunn - LogSumExp_cpp(logunn);
    }

    XI_ijl(j) = exp(tempres);

  }

  return(XI_ijl);
}





// overCAM - specific
arma::mat Update_RHOjk_cpp_overCAM(arma::field<arma::mat> XI_ijl, // collection of  J objectes: nj*L matrices
                                   arma::colvec alpha_star_k,
                                   arma::mat beta_star_lk,
                                   int const L,
                                   int const J,
                                   int const K){

  arma::mat N_jl(J,L);

  for(int j=0; j<J; j++){
    N_jl.row(j) = arma::sum(XI_ijl(j),0);
  }

  arma::mat unn_log_RHO_jk(J,K);
  arma::mat log_RHO_jk(J,K);
  arma::mat RHO_jk(J,K);
  arma::mat E_ln_omega(L,K);

  for(int k=0; k<K; k++){
    E_ln_omega.col(k) =
      E_log_DIR(beta_star_lk.col(k)); // L x K
  }

  arma::colvec ElnPI_k = E_log_DIR(alpha_star_k);

  arma::mat Z = N_jl * E_ln_omega;

  for(int k = 0; k < K; k++){
    unn_log_RHO_jk.col(k) =  Z.col(k) +  ElnPI_k(k);
  }
  for(int j = 0; j < J; j++){
    log_RHO_jk.row(j) = unn_log_RHO_jk.row(j) - LogSumExp_cpp(unn_log_RHO_jk.row(j));
  }

  return(exp(log_RHO_jk));
}


variational_fiSAN <- function(y,
                              group,
                              prior_param = list(m0 = 0,
                                                 tau0 = .01,
                                                 lambda0 = 3,
                                                 gamma0 = 2,
                                                 hyp_alpha1 = 1,
                                                 hyp_alpha2 = 1,
                                                 alpha = NULL,
                                                 b_dirichlet = 1/vi_param$maxL),
                              vi_param = list(maxL = 50,
                                                 maxK = 20,
                                                 epsilon = 1e-6,
                                                 seed = NULL,
                                                 maxSIM = 1e5,
                                                 warmstart = TRUE,
                                                 verbose = FALSE,
                                                 n_runs = 1)){

  ## List completion --------------------------------------------------------

  ### prior_param list
  if(is.null(vi_param$alpha)){
    alpha <- NULL
  }else{
    alpha <- vi_param$alpha
  }

  m0      <- ifelse(is.null(prior_param$m0), 0, prior_param$m0)
  tau0    <- ifelse(is.null(prior_param$tau0), .01, prior_param$tau0)
  lambda0 <- ifelse(is.null(prior_param$lambda0), 3, prior_param$lambda0)
  gamma0  <- ifelse(is.null(prior_param$gamma0), 2, prior_param$gamma0)

  hyp_alpha1 <- ifelse(is.null(prior_param$hyp_alpha1), 1, prior_param$hyp_alpha1)
  hyp_alpha2 <- ifelse(is.null(prior_param$hyp_alpha2), 1, prior_param$hyp_alpha2)


  ### vi_param list
  L <- ifelse(is.null(vi_param$maxL), 50, vi_param$maxL)
  K <- ifelse(is.null(vi_param$maxK), 20, vi_param$maxK)
  epsilon   <- ifelse(is.null(vi_param$epsilon), 1e-6, vi_param$epsilon)
  maxSIM    <- ifelse(is.null(vi_param$maxSIM), 1e5, vi_param$maxSIM)
  warmstart <- ifelse(is.null(vi_param$warmstart), TRUE, vi_param$warmstart)
  verbose   <- ifelse(is.null(vi_param$verbose), FALSE, vi_param$verbose)

  b_dirichlet   <- ifelse(is.null(prior_param$b_dirichlet),  1/L, prior_param$b_dirichlet)

  if(is.null(vi_param$seed)){vi_param$seed <- round(stats::runif(1,1,10000))}
  # random init
  ## Checks -----------------------------------------------------------------

  if(length(y) != length(group)){
    stop("The number of observations and groups must match")
  }

  # Running the function       ----------------------------------------------

  Nj <- tapply(y, group, length)
  J  <- max(group)


  params <- list(y = y,
                group = group,
                Nj = Nj,
                J = J,
                K = K,
                L = L,
                m0 = m0,
                tau0 = tau0,
                lambda0 = lambda0,
                gamma0 = gamma0,
                b_dirichlet = b_dirichlet,
                epsilon =  epsilon,
                seed =  vi_param$seed,
                n_runs = vi_param$n_runs)

  conc_hyperpar <- c(hyp_alpha1,
                     hyp_alpha2)

  conc_par <- c(alpha)


  if(!is.null(conc_par)) {
    if( !(all(conc_par>0)) ){
      stop("Not all the specified concentration parameters are positive.")
    }
    if( length(conc_par) < 1 ){
      stop("Please provide one value for the concentration parameter.")
    }
    # params$alpha1 <- 1 # a tilde
    params$alpha <- conc_par[1] # b tilde
    random_conc <- FALSE
  }else if(!is.null(conc_hyperpar) ) {
    if( !(all(conc_hyperpar>0)) ){
      stop("Not all the specified concentration hyperparameters are positive.")
    }
    if( length(conc_hyperpar) < 2 ){
      stop("Please provide two values for the concentration hyper-parameters.")
    }
    params$hyp_alpha1 <- conc_hyperpar[1]
    params$hyp_alpha2 <- conc_hyperpar[2]
    random_conc <- TRUE
  }else{
    stop("Please provide a consistent specification for the concentration parameters")
  }


  Y_grouped <- list()
  for(j in 1:J){
    Y_grouped[[j]] <- y[group==j] # this is a list, each element is a vector with observations
  }


  set.seed(vi_param$seed)

  # random init
  if(warmstart){
    suppressWarnings(
    ml  <- stats::kmeans(unlist(Y_grouped),centers = L, algorithm="MacQueen",
                         iter.max = 100)$centers
    )
  }else{
    ml  <- stats::runif(L, min(unlist(Y_grouped)),max(unlist(Y_grouped)))
  }


  kl       <- stats::rgamma(L,1,10)
  lambdal  <- stats::rgamma(L,1,1)
  gammal   <- stats::rgamma(L,1,1)



  ###############################################################################################
  XI_ijl <- list()
  for(j in 1:J){
    log.XI_il  <- array(stats::rbeta( Nj[j] * L, 1, 1),dim = c( Nj[j], L))
    Z           <- apply(log.XI_il, c(1), function(x) matrixStats::logSumExp(x))
    XI_ijl[[j]] <- exp(sapply(1:L, function(qq) log.XI_il[,qq]-Z,simplify = "array"))
  }
  ###############################################################################################
  if(K<J){
    log.RHO_jk <- matrix(stats::rbeta(J*K,1,1),J,K)
    Z2         <- apply(log.RHO_jk, c(1), function(x) matrixStats::logSumExp(x))
    RHO_jk     <- exp(sapply(1:K, function(qq) log.RHO_jk[,qq]-Z2,simplify = "matrix"))
  }else if(K==J){
    RHO_jk     <- diag(J)
  }else{
    RHO_jk     <- cbind(diag(J), matrix(0,J,K-J))
  }

  if(random_conc){

    start <- Sys.time()
    results <- main_vb_fiSAN_CP_cpp(
      Y_grouped = Y_grouped,
      XI_ijl = XI_ijl,
      L = L,
      K = K,
      J = J,
      RHO_jk = RHO_jk,
      Nj = Nj,
      m0 = m0,
      k0 = tau0,
      a0 = lambda0,
      b0 = gamma0,
      ml = ml,
      kl = kl,
      al = lambdal,
      bl = gammal,
      conc_hyper = conc_hyperpar,
      beta_bar = b_dirichlet,
      epsilon = epsilon,
      maxSIM = maxSIM,
      verbose = verbose
    )
    end <- Sys.time()

  }else{

    start <- Sys.time()
    results <- main_vb_fiSAN_cpp(
      Y_grouped = Y_grouped,
      XI_ijl = XI_ijl,
      L = L,
      K = K,
      J = J,
      RHO_jk = RHO_jk,
      Nj = Nj,
      m0 = m0,
      k0 = tau0,
      a0 = lambda0,
      b0 = gamma0,
      ml = ml,
      kl = kl,
      al = lambdal,
      bl = gammal,
      a_tilde = 1.0,
      b_tilde = params$alpha,
      beta_bar = b_dirichlet,
      epsilon = epsilon,
      maxSIM = maxSIM,
      verbose = verbose
    )
    end <- Sys.time()

  }
  output <- list("model" = "fiSAN",
                 "params" = params,
                 "sim"= results,
                 "time" = end - start)

  structure(output, class = c("SANvi",class(output)))

}


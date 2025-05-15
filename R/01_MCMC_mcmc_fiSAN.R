#' @importFrom stats cor var dist hclust cutree rgamma

sample_fiSAN <- function(y, group,
                         prior_param = list(
                           m0 = 0,
                           tau0 = 0.1,
                           lambda0 = 3,
                           gamma0 = 2,
                           hyp_alpha1 = 1,
                           hyp_alpha2 = 1,
                           alpha = NULL,
                           b_dirichlet = 1/mcmc_param$maxL),
                         mcmc_param = list(
                           nrep = 1000, burn = 500,
                           maxK = 20, maxL = 50,
                           warmstart = TRUE,
                           nclus_start = NULL,
                           mu_start = NULL,
                           sigma2_start = NULL,
                           M_start = NULL,
                           S_start = NULL,
                           alpha_start = NULL,
                           verbose = TRUE,
                           seed = NULL))
{
  group <- .relabel(group) - 1

  ## List completion --------------------------------------------------------

  ### prior_param list
  prior_param$m0      <- ifelse(is.null(prior_param$m0), 0, prior_param$m0)
  prior_param$tau0    <- ifelse(is.null(prior_param$tau0), .01, prior_param$tau0)
  prior_param$lambda0 <- ifelse(is.null(prior_param$lambda0), 3, prior_param$lambda0)
  prior_param$gamma0  <- ifelse(is.null(prior_param$gamma0), 2, prior_param$gamma0)

  prior_param$hyp_alpha1 <- ifelse(is.null(prior_param$hyp_alpha1), 1, prior_param$hyp_alpha1)
  prior_param$hyp_alpha2 <- ifelse(is.null(prior_param$hyp_alpha2), 1, prior_param$hyp_alpha2)

  ### mcmc_param list
  mcmc_param$nrep <- ifelse(is.null(mcmc_param$nrep), 1000, mcmc_param$nrep)
  mcmc_param$burn <- ifelse(is.null(mcmc_param$burn), 500, mcmc_param$burn)
  mcmc_param$maxL <- ifelse(is.null(mcmc_param$maxL), 50, mcmc_param$maxL)
  mcmc_param$maxK <- ifelse(is.null(mcmc_param$maxK), 20, mcmc_param$maxK)
  mcmc_param$warmstart <- ifelse(is.null(mcmc_param$warmstart), TRUE, mcmc_param$warmstart)
  mcmc_param$verbose <- ifelse(is.null(mcmc_param$verbose), TRUE, mcmc_param$verbose)

  prior_param$b_dirichlet  <- ifelse(is.null(prior_param$b_dirichlet), 1/mcmc_param$maxL, prior_param$b_dirichlet)

  if(is.null(mcmc_param$seed)){mcmc_param$seed <- round(stats::runif(1,1,10000))}
  # random init
  set.seed(mcmc_param$seed)
  ## Checks -----------------------------------------------------------------

  if(length(y) != length(group)){
    stop("The number of observations and groups must match")
  }

  #----------------------------------------------------
  warmstart = mcmc_param$warmstart
  nclus_start = mcmc_param$nclus_start
  mu_start = mcmc_param$mu_start
  sigma2_start = mcmc_param$sigma2_start
  M_start = mcmc_param$M_start
  S_start = mcmc_param$S_start
  verbose = mcmc_param$verbose
  seed = mcmc_param$seed
  #----------------------------------------------------


  params <- list(y = y,
                 group = group+1,
                 Nj = tapply(y,group, length),
                 maxK = mcmc_param$maxK,
                 maxL = mcmc_param$maxL,
                 m0 = prior_param$m0,
                 tau0 = prior_param$tau0,
                 lambda0 = prior_param$lambda0,
                 gamma0 = prior_param$gamma0,
                 b_dirichlet = prior_param$b_dirichlet,
                 seed = seed,
                 nrep = mcmc_param$nrep,
                 burn = mcmc_param$burn)

  if(!is.null(prior_param$alpha)) { params$alpha <- alpha }
  if(is.null(prior_param$alpha)) { params$hyp_alpha1 <- prior_param$hyp_alpha1
                                   params$hyp_alpha2 <- prior_param$hyp_alpha2 }

  if(is.null(S_start)) { S_start <- rep(0,length(unique(group))) }

  # if the initial cluster allocation is passed
  if(!is.null(M_start)) {
    warmstart <- FALSE
    M_start <- .relabel(M_start)

    # and the mean is passed or the variance is passed don't do anything

    # if the mean is not passed
    if(is.null(mu_start)) {
      mu_start <- rep(0,mcmc_param$maxL)
      ncl0 <- length(unique(M_start))
      for(l in unique(M_start)) {
        mu_start[l] <- mean(y[M_start == l])
      }
    }
    # if the variance is not passed
    if(is.null(sigma2_start)) {
      sigma2_start <- rep(0.001,mcmc_param$maxL)
      ncl0 <- length(unique(M_start))
      for(l in unique(M_start)) {
        sigma2_start[l] = var(y[M_start == l])
      }
    }
  } else {
    # if the initial cluster allocation is not passed
    # and you don't want a warmstart
    if(!warmstart){
      M_start <- rep(1, length(y))#sample(0:(maxL-2), length(y), replace = TRUE)
      mu_start <- rep(0, mcmc_param$maxL)
      mu_start[1] <- mean(y)
      sigma2_start <- rep(0.001, mcmc_param$maxL)
      sigma2_start[1] <- var(y)/2
    }

    # if the initial cluster allocation is not passed
    # and you want a warmstart
    if(warmstart){
      mu_start <- rep(0,mcmc_param$maxL)
      sigma2_start <- rep(0.001,mcmc_param$maxL)

      if(is.null(nclus_start)) { nclus_start <- min(c(mcmc_param$maxL, 30))}
      suppressWarnings(
        M_start <- stats::kmeans(y,
                               centers = nclus_start,
                               algorithm="MacQueen",
                               iter.max = 100)$cluster
      )
      nclus_start <- length(unique(M_start))
      mu_start[1:nclus_start] <- sapply(1:nclus_start, function(x) mean(y[M_start == x]))
      sigma2_start[1:nclus_start] <- sapply(1:nclus_start, function(x) var(y[M_start == x]))
      sigma2_start[1:nclus_start][sigma2_start[1:nclus_start]==0] <- 0.001
      sigma2_start[is.na(sigma2_start)] <- 0.001
    }
  }
  M_start <- M_start-1
  sigma2_start[is.na(sigma2_start)] <- 0.001

  fixed_alpha <- F
  if(!is.null(prior_param$alpha) ) {
    fixed_alpha <- T ;
    alpha_start <- prior_param$alpha
  } else { prior_param$alpha <- 1 }


  start <- Sys.time()
  out <- sample_fiSAN_cpp(nrep = mcmc_param$nrep,
                          burn = mcmc_param$burn,
                          y = y, group = group,
                          maxK = mcmc_param$maxK,
                          maxL =  mcmc_param$maxL,
                          m0 = prior_param$m0,
                          tau0 = prior_param$tau0,
                          lambda0 = prior_param$lambda0,
                          gamma0 = prior_param$gamma0,
                          alpha = prior_param$alpha,
                          beta = prior_param$b_dirichlet,
                          hyp_alpha1 = prior_param$hyp_alpha1,
                          hyp_alpha2 = prior_param$hyp_alpha2,
                          fixed_alpha = fixed_alpha,
                          mu_start = mu_start, sigma2_start = sigma2_start,
                          M_start = M_start, S_start = S_start,
                          progressbar = verbose)
  end <- Sys.time()

  warnings <- out$warnings

  out$distr_cluster <- out$distr_cluster + 1
  out$obs_cluster <- out$obs_cluster + 1

  if(length(warnings) == 1) {
    output <- list( "model" = "fiSAN",
                    "params" = params,
                    "sim" = out,
                    "time" = end - start,
                    "warnings" = warnings)
    warning("Increase maxK: all the provided distributional mixture components were used. Check '$warnings' to see when it happened.")
  } else {
    output <- list( "model" = "fiSAN",
                    "params" = params,
                    "sim" = out,
                    "time" = end - start )
  }

  structure(output, class = c("SANmcmc",class(output)))

}

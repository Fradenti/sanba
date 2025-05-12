#' Mean Field Variational Bayes estimation of CAM
#'
#' @noRd
#' @noMd
#'
#' @description \code{variational_CAM} is used to perform posterior inference under the common atoms model (CAM) of Denti et al. (2023) with Gaussian likelihood.
#' The model uses Dirichlet process mixtures (DPM) at both the observational and distributional levels.
#'
#' @usage
#' variational_CAM = function(y,
#'                            group,
#'                            prior_param = list(m0 = 0,
#'                                               tau0 = .01,
#'                                               lambda0 = 3,
#'                                               gamma0 = 2,
#'                                               alpha = NULL,
#'                                               beta = NULL,
#'                                               hyp_alpha1 = 1,
#'                                               hyp_alpha2 = 1,
#'                                               hyp_beta1 = 1,
#'                                               hyp_beta2 = 1),
#'                            vi_param = list(maxL = 30,
#'                                            maxK = 20,
#'                                            epsilon = 1e-6,
#'                                            seed = NULL,
#'                                            maxSIM = 1e5,
#'                                            warmstart = TRUE,
#'                                            verbose = FALSE)
#' )
#'
#' @param y Numerical vector of observations (required).
#' @param group Numerical vector of the same length of \code{y}, indicating the group membership (required).
#' @param maxL,maxK integers, the upper bounds for the observational and distributional clusters to fit, respectively.
#' @param m0,tau0,lambda0,gamma0 Hyperparameters on \eqn{(\mu, \sigma^2) \sim NIG(m_0, \tau_0, \lambda_0,\gamma_0)}.
#' @param conc_hyperpar,conc_par Vectors of values used for the concentration parameters of  of the stick-breaking representation for the distributional and observational DPs, respectively.
#' The following two arguments can be passed. Specifically,
#' \describe{
#'   \item{\code{conc_hyperpar}}{a vector with 4 positive entries: \eqn{(s_1^\alpha,s_2^\alpha,s_1^\beta,s_2^\beta)}.
#'   If a random concentration parameters \eqn{\alpha} and \eqn{\beta} are adopted, the specifications are
#'  \eqn{\alpha \sim Gamma(s_1^\alpha,s_2^\alpha)} and \eqn{\beta \sim Gamma(s_1^\beta,s_2^\beta)}. Default set to unitary vector.}
#'   \item{\code{conc_par}}{a vector with 2 positive entries: \eqn{(\alpha,\beta)}. Default is set to \code{NULL}. If specified, the previous argument is ignored
#'   and the two concentration parameters are assumed fixed and equal to \code{(alpha,beta)}.}
#'  }
#' @param epsilon the tolerance that drives the convergence criterion adopted as stopping rule
#' @param seed random seed to control the initialization.
#' @param maxSIM the maximum number of CAVI iteration to perform.
#' @param warmstart logical, if \code{TRUE}, the observational means of the cluster atoms are initialized with a k-means algorithm.
#' @param verbose logical, if \code{TRUE} the iterations are printed.
#'
#' @details
#'
#' The common atoms mixture model is used to perform inference in nested settings, where the data are organized into \eqn{J} groups.
#' The data should be continuous observations \eqn{(Y_1,\dots,Y_J)}, where each \eqn{Y_j = (y_{1,j},\dots,y_{n_j,j})}
#' contains the \eqn{n_j} observations from group \eqn{j}, for \eqn{j=1,\dots,J}.
#' The function takes as input the data as a numeric vector \code{y} in this concatenated form. Hence \code{y} should be a vector of length
#' \eqn{n_1+\dots+n_J}. The \code{group} parameter is a numeric vector of the same size as \code{y} indicating the group membership for each
#' individual observation.
#' Notice that with this specification the observations in the same group need not be contiguous as long as the correspondence between the variables
#' \code{y} and \code{group} is maintained.
#'
#' \strong{Model}
#'
#' The data are modeled using a Gaussian likelihood, where both the mean and the variance are observational-cluster-specific, i.e.,
#' \deqn{y_{i,j}\mid M_{i,j} = l \sim N(\mu_l,\sigma^2_l)}
#' where \eqn{M_{i,j} \in \{1,2,\dots\}} is the observational cluster indicator of observation \eqn{i} in group \eqn{j}.
#' The prior on the model parameters is a Normal-Inverse-Gamma distribution \eqn{(\mu_l,\sigma^2_l)\sim NIG (m_0,\tau_0,\lambda_0,\gamma_0)},
#' i.e., \eqn{\mu_l\mid\sigma^2_l \sim N(m_0, \sigma^2_l / \tau_0)}, \eqn{1/\sigma^2_l \sim Gamma(\lambda_0, \gamma_0)} (shape, rate).
#'
#' \strong{Clustering}
#'
#' The model performs a clustering of both observations and groups.
#' The clustering of groups (distributional clustering) is provided by the allocation variables \eqn{S_j \in \{1,2,\dots\}}, with
#' \deqn{Pr(S_j = k \mid \dots ) = \pi_k  \qquad \text{for } \: k = 1,2,\dots}
#' The distribution of the probabilities is \eqn{ \{\pi_k\}_{k=1}^{\infty} \sim GEM(\alpha)},
#' where GEM is the Griffiths-Engen-McCloskey distribution of parameter \eqn{\alpha},
#' which characterizes the stick-breaking construction of the DP (Sethuraman, 1994).
#'
#' The clustering of observations (observational clustering) is provided by the allocation variables \eqn{M_{i,j} \in \{1,2,\dots\}}, with
#' \deqn{ Pr(M_{i,j} = l \mid S_j = k, \dots ) = \omega_{l,k} \qquad \text{for } \: k = 1,2,\dots \, ; \: l = 1,2,\dots }
#' The distribution of the probabilities is \eqn{ \{\omega_{l,k}\}_{l=1}^{\infty} \sim GEM(\beta)} for all \eqn{k = 1,2,\dots}
#'
#'
#'
#' @return \code{variational_CAM} returns a list of class \code{SANvi} containing four objects:
#' \itemize{
#'   \item \code{model}: name of the fitted model.
#'   \item \code{params}: list containing the data and the parameters used in the simulation. Details below.
#'   \item \code{sim}: list containing the simulated values (optimized variational parameters). Details below.
#'   \item \code{time}: total computation time.
#' }
#'
#'
#' \strong{Data and parameters}:
#' \code{params} is a list with the following components:
#' \describe{
#' \item{\code{y, group, Nj, J}}{Data, group labels, group frequencies, and number of groups.}
#' \item{\code{K, L}}{Number of fitted distributional and observational clusters.}
#' \item{\code{m0, tau0, lambda0, gamma0}}{Model hyperparameters.}
#' \item{\code{epsilon, seed}}{The threshold controlling the convergence criterion and the random seed adopted to replicate the run.}
#' \item{(\code{hyp_alpha1,hyp_alpha2}) or \code{alpha}}{Hyperparameters on \eqn{\alpha} (if \eqn{\alpha} random);
#'   or provided value for \eqn{\alpha} (if fixed).}
#' \item{(\code{hyp_beta1,hyp_beta2}) or \code{beta}}{Hyperparameters on \eqn{\beta} (if \eqn{\beta} random);
#'   or provided value for \eqn{\beta} (if fixed).}
#' }
#'
#'
#' \strong{Simulated values}:
#' \code{sim} is a list with the following components:
#' \describe{
#' \item{\code{theta_l}}{Matrix of size (L,4).
#'    Each row is a posterior variational estimate of the four normal-inverse gamma hyperparameters.}
#' \item{\code{Elbo_val}}{Vector containing the values of the ELBO.}
#' \item{\code{XI}}{A list of length J. Each element is a matrix of size (N, L)
#'    posterior variational probability of assignment of assignment of the i-th observation in the j-th group to the l-th OC,
#'    i.e., \eqn{\hat{\xi}_{i,j,l} = \hat{\mathbb{Q}}(M_{i,j}=l)}.}
#' \item{\code{RHO}}{Matrix of size (J, K).
#'    Each row is a posterior variational probability of assignment of the j-th group to the k-th DC, i.e., \eqn{\hat{\rho}_{j,k} = \hat{\mathbb{Q}}(S_j=k)}. }
#' \item{\code{a_tilde_k,b_tilde_k}}{Vector of updated variational parameters of the Beta distributions governing the distributional stick-breaking process.}
#' \item{\code{a_tilde_lk,b_tilde_lk}}{Matrix of updated variational parameters of the Beta distributions governing the observational stick-breaking
#'                        process (arranged by column).}
#' \item{\code{conc_hyper}}{If the concentration parameters are chosen to be random, these object contain a vector with the four updated hyperparameters.}
#' \item{\code{alpha},\code{beta}}{If the concentration parameters are chosen to be fixed, these objects contain the passed values.}
#'}
#'
#'
#'
#' @references Denti, F., Camerlenghi, F., Guindani, M., and Mira, A. (2023). A Common Atoms Model for the Bayesian Nonparametric Analysis of Nested Data.
#' \emph{Journal of the American Statistical Association}, 118(541), 405-416. DOI: 10.1080/01621459.2021.1933499
#'
#' Sethuraman, A.J. (1994). A Constructive Definition of Dirichlet Priors, \emph{Statistica Sinica}, 4, 639–650.
#'
#'
#'
#'
variational_CAM = function(y,
                  group,
                  prior_param = list(m0 = 0,
                                     tau0 = .01,
                                     lambda0 = 3,
                                     gamma0 = 2,
                                     alpha = NULL,
                                     beta = NULL,
                                     hyp_alpha1 = 1,
                                     hyp_alpha2 = 1,
                                     hyp_beta1 = 1,
                                     hyp_beta2 = 1),
                  vi_param = list(maxL = 30,
                                     maxK = 20,
                                     epsilon = 1e-6,
                                     seed = NULL,
                                     maxSIM = 1e5,
                                     warmstart = TRUE,
                                     verbose = FALSE,
                                     n_runs = 1)
                  )
{

  # Checks and list completion ----------------------------------------------


  ## List completion --------------------------------------------------------

  ### prior_param list
  if(is.null(vi_param$alpha)){
    alpha <- NULL
  }else{
    alpha <- vi_param$alpha
  }
  if(is.null(vi_param$beta)){
    beta <- NULL
  }else{
    beta <- vi_param$beta
  }

  m0      <- ifelse(is.null(vi_param$m0),   0, vi_param$m0)
  tau0    <- ifelse(is.null(vi_param$tau0), .01, vi_param$tau0)
  lambda0 <- ifelse(is.null(vi_param$lambda0), 3, vi_param$lambda0)
  gamma0  <- ifelse(is.null(vi_param$gamma0), 2, vi_param$gamma0)

  hyp_alpha1 <- ifelse(is.null(vi_param$hyp_alpha1), 1, vi_param$hyp_alpha1)
  hyp_alpha2 <- ifelse(is.null(vi_param$hyp_alpha2), 1, vi_param$hyp_alpha2)
  hyp_beta1  <- ifelse(is.null(vi_param$hyp_beta1),  1, vi_param$hyp_beta1)
  hyp_beta2  <- ifelse(is.null(vi_param$hyp_beta2),  1, vi_param$hyp_beta2)


  ### vi_paramm list
  L <- ifelse(is.null(vi_param$maxL), 30, vi_param$maxL)
  K <- ifelse(is.null(vi_param$maxK), 20, vi_param$maxK)
  epsilon   <- ifelse(is.null(vi_param$epsilon), 1e-6,   vi_param$epsilon)
  maxSIM    <- ifelse(is.null(vi_param$maxSIM), 1e5,     vi_param$maxSIM)
  warmstart <- ifelse(is.null(vi_param$warmstart), TRUE, vi_param$warmstart)
  verbose   <- ifelse(is.null(vi_param$verbose), FALSE,  vi_param$verbose)

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
                seed = vi_param$seed,
                epsilon = epsilon,
                n_runs = vi_param$n_runs
                )

  conc_hyperpar <- c(hyp_alpha1,
                     hyp_alpha2,
                     hyp_beta1,
                     hyp_beta2)

  conc_par <- c(alpha, beta)

  if(!is.null(conc_par)) {
    if( !(all(conc_par>0))){
      stop("Not all the specified concentration parameters are positive.")
    }
    if( length(conc_par) < 2 ){
      stop("Please provide two values for the concentration parameters.")
    }
    # params$alpha1 <- 1 # a tilde
    params$alpha <- conc_par[1] # b tilde
    # params$beta1  <- 1  # a bar
    params$beta  <- conc_par[2]  # b bar
    random_conc <- FALSE
  }else if(!is.null(conc_hyperpar)) {
    if( !(all(conc_hyperpar>0) )){
      stop("Not all the specified concentration hyperparameters are positive.")
    }
    if( length(conc_hyperpar) < 4 ){
      stop("Please provide four values for the concentration hyper-parameters.")
    }

    params$hyp_alpha1 <- conc_hyperpar[1]
    params$hyp_alpha2 <- conc_hyperpar[2]
    params$hyp_beta1  <- conc_hyperpar[3]
    params$hyp_beta2  <- conc_hyperpar[4]
    random_conc <- TRUE
  }else{
    stop("Please provide a consistent specification for the concentration parameters")
  }


  Y_grouped <- list()
  for(j in 1:J){
    Y_grouped[[j]] <- y[group==j] # this is a list, each element is a vector with observations
  }


  set.seed(vi_param$seed)


  if(warmstart){
    suppressWarnings(
      ml  <- stats::kmeans(unlist(Y_grouped),
                         centers = L,
                         algorithm="MacQueen",
                         iter.max = 100)$centers
    )
  }else{
    ml  <- stats::runif(L, min(unlist(Y_grouped)),max(unlist(Y_grouped)))
  }

  kl       <- stats::rgamma(L,1,10)
  lambdal  <- stats::rgamma(L,1,1)
  gammal   <- stats::rgamma(L,1,1)


  ###############################################################################################
  XI_ijl = list()
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
    results <- main_vb_cam_CP_cpp(
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
     epsilon = epsilon,
     maxSIM =  maxSIM,
     verbose = verbose
    )
    end <- Sys.time()

  }else{

    start <- Sys.time()
    results <- main_vb_cam_cpp(
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
      a_tilde = 1,
      b_tilde = params$alpha,
      a_bar = 1,
      b_bar = params$beta,
      epsilon = epsilon,
      maxSIM =  maxSIM,
      verbose = verbose
    )
    end <- Sys.time()
  }

  output <- list("model" = "CAM",
                 "params" = params,
                 "sim"= results,
                 "time" = end - start)

  structure(output, class = c("SANvi",class(output)))

}



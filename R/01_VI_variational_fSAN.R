#' Mean Field Variational Bayes estimation of fSAN
#'
#' @noRd
#' @noMd
#'
#' @description \code{variational_fSAN} is used to perform posterior inference under the finite shared atoms nested (fSAN) model with Gaussian likelihood (originally proposed in D'Angelo et al., 2023).
#' The model uses finite Dirichlet mixtures for both the distributional and observational levels of the model.
#'
#'
#' @usage
#' variational_fSAN(y,
#'                  group,
#'                  prior_param = list(m0 = 0,
#'                                     tau0 = .01,
#'                                     lambda0 = 3,
#'                                     gamma0 = 2,
#'                                     a_dirichlet = .005,
#'                                     b_dirichlet = .005),
#'                  vi_param = list(maxL = 30,
#'                                  maxK = 20,
#'                                  epsilon = 1e-6,
#'                                  seed = NULL,
#'                                  maxSIM = 1e5,
#'                                  warmstart = TRUE,
#'                                  verbose = FALSE))
#'
#' @param y Numerical vector of observations (required).
#' @param group Numerical vector of the same length of \code{y}, indicating the group membership (required).
#' @param maxL,maxK integers, the upper bounds for the observational and distributional clusters to fit, respectively
#' @param m0,tau0,lambda0,gamma0 Hyperparameters on \eqn{(\mu, \sigma^2) \sim NIG(m_0, \tau_0, \lambda_0,\gamma_0)}.
#' @param a_dirichlet the hyperparameter of the symmetric distributional Dirichlet distribution.
#' @param b_dirichlet the hyperparameter of the symmetric observational Dirichlet distribution.
#' @param epsilon the tolerance that drives the convergence criterion adopted as stopping rule
#' @param seed random seed to control the initialization.
#' @param maxSIM the maximum number of CAVI iteration to perform.
#' @param warmstart logical, if \code{TRUE}, the observational means of the cluster atoms are initialized with a k-means algorithm.
#' @param verbose logical, if \code{TRUE} the iterations are printed.
#'
#' @details
#' \strong{Data structure}
#'
#' The finite common atoms mixture model is used to perform inference in nested settings, where the data are organized into \eqn{J} groups.
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
#' where \eqn{M_{i,j} \in \{1,\dots,L \}} is the observational cluster indicator of observation \eqn{i} in group \eqn{j}.
#' The prior on the model parameters is a Normal-Inverse-Gamma distribution \eqn{(\mu_l,\sigma^2_l)\sim NIG (m_0,\tau_0,\lambda_0,\gamma_0)},
#' i.e., \eqn{\mu_l\mid\sigma^2_l \sim N(m_0, \sigma^2_l / \tau_0)}, \eqn{1/\sigma^2_l \sim Gamma(\lambda_0, \gamma_0)} (shape, rate).
#'
#' \strong{Clustering}
#'
#' The model performs a clustering of both observations and groups.
#' The clustering of groups (distributional clustering) is provided by the allocation variables \eqn{S_j \in \{1,\dots,K\}}, with
#' \deqn{Pr(S_j = k \mid \dots ) = \pi_k  \qquad \text{for } \: k = 1,\dots,K.}
#' The distribution of the probabilities is \eqn{(\pi_1,\dots,\pi_{K})\sim Dirichlet_K(\alpha/K,\dots,\alpha/K)}.
#' Here, the dimension \eqn{K} is fixed.
#'
#' The clustering of observations (observational clustering) is provided by the allocation variables \eqn{M_{i,j} \in \{1,\dots,L\}}, with
#' \deqn{ Pr(M_{i,j} = l \mid S_j = k, \dots ) = \omega_{l,k} \qquad \text{for } \: k = 1,\dots,K \, ; \: l = 1,\dots,L. }
#' The distribution of the probabilities is \eqn{(\omega_{1,k},\dots,\omega_{L,k})\sim Dirichlet_L(\beta/L,\dots,\beta/L)} for all \eqn{k = 1,\dots,K}.
#' Here, the dimension \eqn{L} is fixed.
#'
#' @return \code{variational_fSAN} returns a list of class \code{SANvi} containing four objects:
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
#' \item{\code{a_dirichlet, b_dirichlet}}{the hyperparameters governing all the finite Dirichlet distributions at the distributional and observational level.}
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
#' \item{\code{a_dirichlet_k}}{Vector of updated variational parameters of the Dirichlet distributions governing the distributional clustering.}
#' \item{\code{b_dirichlet_lk}}{Matrix of updated variational parameters of the Dirichlet distributions governing the observational clustering (arranged by column).}
#'}
#'
#'
#' @references D’Angelo, L., Canale, A., Yu, Z., and Guindani, M. (2023).
#' Bayesian nonparametric analysis for the detection of spikes in noisy calcium imaging data. \emph{Biometrics}, 79(2), 1370–1382. DOI: 10.1111/biom.13626
#'
#'
variational_fSAN <- function(y,
                    group,
                    prior_param = list(m0 = 0,
                                       tau0 = .01,
                                       lambda0 = 3,
                                       gamma0 = 2,
                                       a_dirichlet = .005,
                                       b_dirichlet = .005),
                    vi_param = list(maxL = 30,
                                       maxK = 20,
                                       epsilon = 1e-6,
                                       seed = NULL,
                                       n_runs = 1,
                                       maxSIM = 1e5,
                                       warmstart = TRUE,
                                       verbose = FALSE,
                                       n_runs = 1)
)
{
  ## List completion --------------------------------------------------------

  ### prior_param list
  m0      <- ifelse(is.null(prior_param$m0), 0, prior_param$m0)
  tau0    <- ifelse(is.null(prior_param$tau0), .01, prior_param$tau0)
  lambda0 <- ifelse(is.null(prior_param$lambda0), 3, prior_param$lambda0)
  gamma0  <- ifelse(is.null(prior_param$gamma0), 2, prior_param$gamma0)

  a_dirichlet <- ifelse(is.null(prior_param$a_dirichlet), .005, prior_param$a_dirichlet)
  b_dirichlet   <- ifelse(is.null(prior_param$b_dirichlet), .005, prior_param$b_dirichlet)


  ### vi_param list
  L <- ifelse(is.null(vi_param$maxL), 30, vi_param$maxL)
  K <- ifelse(is.null(vi_param$maxK), 20, vi_param$maxK)
  epsilon   <- ifelse(is.null(vi_param$epsilon), 1e-6, vi_param$epsilon)
  maxSIM    <- ifelse(is.null(vi_param$maxSIM), 1e5, vi_param$maxSIM)
  warmstart <- ifelse(is.null(vi_param$warmstart), TRUE, vi_param$warmstart)
  verbose   <- ifelse(is.null(vi_param$verbose), FALSE, vi_param$verbose)

  if(is.null(vi_param$seed)){vi_param$seed <- round(stats::runif(1,1,10000))}

  ## Checks -----------------------------------------------------------------

  if(length(y) != length(group)){
    stop("The number of observations and groups must match")
  }

  # Running the function       ----------------------------------------------

  Nj <- tapply(y,group, length)
  J  <- max(group)

  if(is.null( vi_param$seed)){ vi_param$seed <- round(stats::runif(1,1,10000))}
  # random init


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
                a_dirichlet = a_dirichlet,
                b_dirichlet = b_dirichlet,
                epsilon =  vi_param$epsilon,
                seed =  vi_param$seed,
                n_runs = vi_param$n_runs)

  Y_grouped <- list()
  for(j in 1:J){
    Y_grouped[[j]] <- y[group==j] # this is a list, each element is a vector with observations
  }

  # random init
  set.seed(vi_param$seed)

  if(warmstart){
    suppressWarnings(
      ml  <- stats::kmeans(unlist(Y_grouped),centers = L, algorithm="MacQueen",
                         iter.max = 100)$centers
    )
  }else{
    ml  <- stats::runif(L, min(unlist(Y_grouped)),max(unlist(Y_grouped)))
  }
  kl  <- stats::rgamma(L,1,10)
  lambdal  <- stats::rgamma(L,1,1)
  gammal  <- stats::rgamma(L,1,1)


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



  start <- Sys.time()
  results <- main_vb_fSAN_cpp(
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
    alpha_bar =  a_dirichlet,
    beta_bar =   b_dirichlet,
    epsilon = epsilon,
    maxSIM =  maxSIM,
    verbose = verbose
  )
  end <- Sys.time()

  output <- list("model" = "fSAN",
                 "params" = params,
                 "sim"= results,
                 "time" = end - start)

  structure(output, class = c("SANvi",class(output)))

}




#' Fit the Finite Shared Atoms Mixture Model
#'
#' @description \code{fit_fSAN} fits the finite shared atoms nested (fSAN) mixture model with Gaussian kernels and normal-inverse gamma priors on the unknown means and variances.
#' The function returns an object of class \code{SANmcmc} or \code{SANvi} depending on the chosen computational approach (MCMC or VI).
#'
#' @usage
#' fit_fSAN(y, group, est_method = c("VI", "MCMC"),
#'          prior_param = list(),
#'          vi_param = list(),
#'          mcmc_param = list())
#'
#' @param y Numerical vector of observations (required).
#' @param group Numerical vector of the same length of \code{y}, indicating the group membership (required).
#' @param est_method Character, specifying the preferred estimation method. It can be either \code{"VI"} or \code{"MCMC"}.
#' @param prior_param A list containing
#' \describe{
#'    \item{\code{m0, tau0, lambda0, gamma0}}{Hyperparameters on \eqn{(\mu, \sigma^2) \sim NIG(m_0, \tau_0, \lambda_0,\gamma_0)}. The default is (0, 0.01, 3, 2).}
#'    \item{\code{a_dirichlet}}{The hyperparameter of the symmetric distributional Dirichlet distribution. The default is 1/\code{maxK}.}
#'    \item{\code{b_dirichlet}}{The hyperparameter of the symmetric observational Dirichlet distribution. The default is 1/\code{maxL}.}
#' }
#'
#' @param vi_param A list of variational inference-specific settings, containing
#'  \describe{
#'    \item{\code{maxL, maxK}}{Integers, the upper bounds for the observational and distributional clusters to fit, respectively. The default is (50, 20).}
#'    \item{\code{epsilon}}{The threshold controlling the convergence criterion.}
#'    \item{\code{n_runs}}{Number of starting points considered for the estimation.}
#'    \item{\code{seed}}{Random seed to control the initialization.}
#'    \item{\code{maxSIM}}{The maximum number of CAVI iteration to perform.}
#'    \item{\code{warmstart}}{Logical, if \code{TRUE}, the observational means of the cluster atoms are initialized with a k-means algorithm.}
#'    \item{\code{verbose}}{Logical, if \code{TRUE} the iterations are printed.}
#' }
#' @param mcmc_param A list of MCMC inference-specific settings, containing
#'  \describe{
#'    \item{\code{nrep, burn}}{Integers, the number of total MCMC iterations, and the number of discarded iterations, respectively.}
#'    \item{\code{maxL, maxK}}{Integers, the upper bounds for the observational and distributional clusters to fit, respectively. The default is (50, 20).}
#'    \item{\code{seed}}{Random seed to control the initialization.}
#'    \item{\code{warmstart}}{Logical, if \code{TRUE}, the observational means of the cluster atoms are initialized with a k-means algorithm. If \code{FALSE}, the starting points can be passed through the parameters \code{ nclus_start, mu_start, sigma2_start, M_start, S_start} }
#'    \item{\code{verbose}}{Logical, if \code{TRUE} the iterations are printed.}
#' }
#'
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
#' The prior on the model parameters is a normal-inverse gamma distribution \eqn{(\mu_l,\sigma^2_l)\sim NIG (m_0,\tau_0,\lambda_0,\gamma_0)},
#' i.e., \eqn{\mu_l\mid\sigma^2_l \sim N(m_0, \sigma^2_l / \tau_0)}, \eqn{1/\sigma^2_l \sim Gamma(\lambda_0, \gamma_0)} (shape, rate).
#'
#' \strong{Clustering}
#'
#' The model performs a clustering of both observations and groups.
#' The clustering of groups (distributional clustering) is provided by the allocation variables \eqn{S_j \in \{1,\dots,K\}}, with
#' \deqn{Pr(S_j = k \mid \dots ) = \pi_k  \qquad \text{for } \: k = 1,\dots,K.}
#' The distribution of the probabilities is \eqn{(\pi_1,\dots,\pi_{K})\sim Dirichlet_K(a,\dots,a)}.
#' Here, the dimension \eqn{K} is fixed.
#'
#' The clustering of observations (observational clustering) is provided by the allocation variables \eqn{M_{i,j} \in \{1,\dots,L\}}, with
#' \deqn{ Pr(M_{i,j} = l \mid S_j = k, \dots ) = \omega_{l,k} \qquad \text{for } \: k = 1,\dots,K \, ; \: l = 1,\dots,L. }
#' The distribution of the probabilities is \eqn{(\omega_{1,k},\dots,\omega_{L,k})\sim Dirichlet_L(b,\dots,b)} for all \eqn{k = 1,\dots,K}.
#' Here, the dimension \eqn{L} is fixed.
#'
#' @return \code{fit_fSAN} returns a list of class \code{SANvi}, if \code{est_method = "VI"}, or \code{SANmcmc}, if \code{est_method = "MCMC"}. The list contains the following elements:
#' \describe{
#'   \item{ \code{model}}{Name of the fitted model.}
#'   \item{ \code{params}}{List containing the data and the parameters used in the simulation. Details below.}
#'   \item{ \code{sim}}{List containing the optimized variational parameters or the simulated values. Details below.}
#'   \item{ \code{time}}{Total computation time.}
#' }
#'
#'
#' \strong{Data and parameters}:
#' \code{params} is a list with the following components:
#' \itemize{
#' \item \code{y, group, Nj, J}: Data, group labels, group frequencies, and number of groups.
#' \item \code{K, L}: Number of distributional and observational mixture components.
#' \item \code{m0, tau0, lambda0, gamma0}: Model hyperparameters.
#' \item \code{a_dirichlet}: Provided value for \eqn{a}.
#' \item \code{b_dirichlet}: Provided value for \eqn{b}.
#' \item \code{seed}: The random seed adopted to replicate the run.
#' \item \code{epsilon, n_runs}: The threshold controlling the convergence criterion and the number of starting points considered
#' \item \code{nrep, burnin}: If \code{est_method = "MCMC"}, the number of total MCMC iterations, and the number of discarded ones.
#' }
#'
#' \strong{Simulated values}: depending on the algorithm, it returns a list with the optimized variational parameters or a list with the chains of the simulated values.
#'
#' \strong{Variational inference}: \code{sim} is a list with the following components:
#' \itemize{
#' \item \code{theta_l}: Matrix of size (\code{maxL}, 4).
#'    Each row is a posterior variational estimate of the four normal-inverse gamma hyperparameters.
#' \item \code{XI} : A list of length J. Each element is a matrix of size (\code{Nj}, \code{maxL})
#'    posterior variational probability of assignment of assignment of the i-th observation in the j-th group to the l-th OC,
#'    i.e., \eqn{\hat{\xi}_{i,j,l} = \hat{\mathbb{Q}}(M_{i,j}=l)}.
#' \item \code{RHO}: Matrix of size (J, \code{maxK}).
#'    Each row is a posterior variational probability of assignment of the j-th group to the k-th DC, i.e., \eqn{\hat{\rho}_{j,k} = \hat{\mathbb{Q}}(S_j=k)}.
#' \item \code{a_dirichlet_k}: Vector of updated variational parameters of the Dirichlet distribution governing the distributional clustering.
#' \item \code{b_dirichlet_lk}: Matrix of updated variational parameters of the Dirichlet distributions governing the observational clustering (arranged by column).
#' \item \code{Elbo_val}: Vector containing the values of the ELBO.
#'}
#'
#' \strong{MCMC inference}: \code{sim} is a list with the following components:
#' \itemize{
#' \item \code{mu}: Matrix of size (\code{nrep}, \code{maxL}).
#'    Each row is a posterior sample of the mean parameter of each observational cluster \eqn{(\mu_1,\dots\mu_L)}.
#' \item \code{sigma2}: Matrix of size (\code{nrep}, \code{maxL}).
#'     Each row is a posterior sample of the variance parameter of each observational cluster \eqn{(\sigma^2_1,\dots\sigma^2_L)}.
#' \item \code{obs_cluster}: Matrix of size (\code{nrep}, n), with n = \code{length(y)}.
#'    Each row is a posterior sample of the observational cluster allocation variables \eqn{(M_{1,1},\dots,M_{n_J,J})}.
#' \item \code{distr_cluster}: Matrix of size (\code{nrep}, J), with J = \code{length(unique(group))}.
#'    Each row is a posterior sample of the distributional cluster allocation variables \eqn{(S_1,\dots,S_J)}.
#' \item \code{pi}: Matrix of size (\code{nrep}, \code{maxK}).
#'    Each row is a posterior sample of the distributional cluster probabilities \eqn{(\pi_1,\dots,\pi_{\code{maxK}})}.
#' \item \code{omega}: 3-d array of size (\code{maxL}, \code{maxK}, \code{nrep}).
#'    Each slice is a posterior sample of the observational cluster probabilities.
#'    In each slice, each column \eqn{k} is a vector (of length \code{maxL}) observational cluster probabilities
#'    \eqn{(\omega_{1,k},\dots,\omega_{\code{maxL},k})} for distributional cluster \eqn{k}.
#'}
#'
#'
#'
#'@export
#'
#'
#'@examples
#' set.seed(123)
#' y <- c(rnorm(60), rnorm(40, 5))
#' g <- rep(1:2, rep(50, 2))
#' plot(density(y[g==1]), xlim = c(-5,10), main = "Group-specific density")
#' lines(density(y[g==2]), col = 2)
#'
#' out_vi <- fit_fSAN(y, group = g, est_method = "VI", vi_param = list(n_runs = 1))
#' out_vi
#'
#' out_mcmc <- fit_fSAN(y = y, group = g, est_method = "MCMC",
#'                       mcmc_param = list(nrep = 100, burn= 50))
#' out_mcmc
#'
fit_fSAN <- function(y,
                     group,
                     est_method = c("VI", "MCMC"),
                     prior_param = list(),
                     vi_param = list(),
                     mcmc_param = list()
                     ){

  # Checks and list completion ----------------------------------------------
  est_method <- match.arg(arg = est_method)


  if(est_method == "VI"){
    vi_param$n_runs    <- ifelse(is.null(vi_param$n_runs), 1, vi_param$n_runs)

    if(vi_param$n_runs == 1){
    est_model <- variational_fSAN(y,
                                  group,
                                  prior_param = prior_param,
                                  vi_param = vi_param)
    }else{
      print_run_progress    <- ifelse(is.null(vi_param$print_run_progress), FALSE, vi_param$print_run_progress)

      list_est <- list()
      elbos    <- list()
      if(is.null( vi_param$seed)){
        ROOT <- round(stats::runif(1,1,10000))
      }else{
        ROOT <- vi_param$seed
      }

      for(l in 1:vi_param$n_runs){

        vi_param$seed <- ROOT*l

        proposed_fit <- variational_fSAN(y,
                                    group,
                                    prior_param = prior_param,
                                    vi_param = vi_param)
        elbos[[l]] <- proposed_fit$sim$Elbo_val
        if( l == 1){
          est_model <- proposed_fit
          max_elbo_observed <- max(elbos[[1]])
        } else if(l > 2) {
          if( max_elbo_observed < max(elbos[[l]]) ){
            est_model <- proposed_fit
            max_elbo_observed <- max(elbos[[l]])
          }
        }

        if(print_run_progress){
          cat(paste0("Performing run number ", l, " out of ", vi_param$n_runs, "\n"))
        }
      }

      est_model$all_elbos <- elbos
    }

    }else{
    est_model <- sample_fSAN(y,
                             group,
                             prior_param = prior_param,
                             mcmc_param = mcmc_param)
  }





  return(est_model)
}

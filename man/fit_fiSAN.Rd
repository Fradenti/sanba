% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_fiSAN.R
\name{fit_fiSAN}
\alias{fit_fiSAN}
\title{Fit the Finite-Infinite Shared Atoms Mixture Model}
\usage{
fit_fiSAN(y, group, est_method = c("VI", "MCMC"),
         prior_param = list(),
         vi_param = list(),
         mcmc_param = list())
}
\arguments{
\item{y}{Numerical vector of observations (required).}

\item{group}{Numerical vector of the same length of \code{y}, indicating the group membership (required).}

\item{est_method}{Character, specifying the preferred estimation method. It can be either \code{"VI"} or \code{"MCMC"}.}

\item{prior_param}{A list containing:
\describe{
  \item{\code{m0, tau0, lambda0, gamma0}}{Hyperparameters on \eqn{(\mu, \sigma^2) \sim NIG(m_0, \tau_0, \lambda_0,\gamma_0)}. The default is (0, 0.01, 3, 2).}
  \item{\code{hyp_alpha1, hyp_alpha2}}{If a random \eqn{\alpha} is used, (\code{hyp_alpha1}, \code{hyp_alpha2}) specify the hyperparameters. The default is (1,1). The prior is \eqn{\alpha} ~ Gamma(\code{hyp_alpha1}, \code{hyp_alpha2}).}
  \item{\code{alpha}}{Distributional DP parameter if fixed (optional). The distribution is \eqn{\pi\sim \text{GEM} (\alpha)}.}
  \item{\code{b_dirichlet}}{The hyperparameter of the symmetric observational Dirichlet distribution. The default is 1/\code{maxL}.}
}}

\item{vi_param}{A list of variational inference-specific settings containing:
\describe{
  \item{\code{maxL, maxK}}{Integers, the upper bounds for the observational and distributional clusters to fit, respectively. The default is (50, 20).}
   \item{\code{epsilon}}{The threshold controlling the convergence criterion.}
   \item{\code{n_runs}}{Number of starting points considered for the estimation.}
  \item{\code{seed}}{Random seed to control the initialization.}
  \item{\code{maxSIM}}{The maximum number of CAVI iterations to perform.}
  \item{\code{warmstart}}{Logical, if \code{TRUE}, the observational means of the cluster atoms are initialized with a k-means algorithm.}
  \item{\code{verbose}}{Logical, if \code{TRUE} the iterations are printed.}
}}

\item{mcmc_param}{A list of MCMC inference-specific settings containing:
\describe{
  \item{\code{nrep, burn}}{Integers, the number of total MCMC iterations, and the number of discarded iterations, respectively.}
  \item{\code{maxL, maxK}}{Integers, the upper bounds for the observational and distributional clusters to fit, respectively. The default is (50, 20).}
  \item{\code{seed}}{Random seed to control the initialization.}
  \item{\code{warmstart}}{Logical, if \code{TRUE}, the observational means of the cluster atoms are initialized with a k-means algorithm. If \code{FALSE}, the starting points can be passed through the parameters \code{nclus_start, mu_start, sigma2_start, M_start, S_start, alpha_start}.}
  \item{\code{verbose}}{Logical, if \code{TRUE} the iterations are printed.}
}}
}
\value{
\code{fit_fiSAN} returns a list of class \code{SANvi}, if \code{est_method = "VI"}, or \code{SANmcmc}, if \code{est_method = "MCMC"}. The list contains the following elements:
\describe{
  \item{\code{model}}{Name of the fitted model.}
  \item{\code{params}}{List containing the data and the parameters used in the simulation. Details below.}
  \item{\code{sim}}{List containing the optimized variational parameters or the simulated values. Details below.}
  \item{\code{time}}{Total computation time.}
}

\strong{Data and parameters}:
\code{params} is a list with the following components:
\itemize{
  \item \code{y, group, Nj, J}: Data, group labels, group frequencies, and number of groups.
  \item \code{K, L}: Number of distributional and observational mixture components.
  \item \code{m0, tau0, lambda0, gamma0}: Model hyperparameters.
  \item (\code{hyp_alpha1, hyp_alpha2}) or \code{alpha}: Hyperparameters on \eqn{\alpha} (if \eqn{\alpha} random);
  or provided value for \eqn{\alpha} (if fixed).
  \item \code{b_dirichlet}: Provided value for \eqn{b}.
  \item \code{seed}: The random seed adopted to replicate the run.
  \item \code{epsilon, n_runs}: If \code{est_method = "VI"}, the threshold controlling the convergence criterion and the number of
  starting points considered.
  \item \code{nrep, burnin}: If \code{est_method = "MCMC"}, the number of total MCMC iterations, and the number of discarded ones.
}

\strong{Simulated values}: Depending on the algorithm, it returns a list with the optimized variational parameters or a list with the chains of the simulated values.

\strong{Variational inference}: \code{sim} is a list with the following components:
\itemize{
  \item \code{theta_l}: Matrix of size (\code{maxL}, 4). Each row is a posterior variational estimate of the four normal-inverse gamma hyperparameters.
  \item \code{XI}: A list of length J. Each element is a matrix of size (\code{Nj}, \code{maxL}), the posterior variational assignment probabilities \eqn{\hat{\mathbb{Q}}(M_{i,j}=l)}.
  \item \code{RHO}: Matrix of size (J, \code{maxK}), with the posterior variational assignment probabilities \eqn{\hat{\mathbb{Q}}(S_j=k)}.
  \item \code{a_tilde_k, b_tilde_k}: Vector of updated variational parameters of the beta distributions governing the distributional stick-breaking process.
  \item \code{conc_hyper}: If the concentration parameter is random, this contains its updated hyperparameters.
  \item \code{b_dirichlet_lk}: Matrix of updated variational parameters of the Dirichlet distributions governing observational clustering.
  \item \code{Elbo_val}: Vector containing the values of the ELBO.
}

\strong{MCMC inference}: \code{sim} is a list with the following components:
\itemize{
  \item \code{mu}: Matrix of size (\code{nrep}, \code{maxL}) with samples of the observational cluster means.
  \item \code{sigma2}: Matrix of size (\code{nrep}, \code{maxL}) with samples of the observational cluster variances.
  \item \code{obs_cluster}: Matrix of size (\code{nrep}, n) with posterior samples of the observational cluster allocations.
  \item \code{distr_cluster}: Matrix of size (\code{nrep}, J) with posterior samples of the distributional cluster allocations.
  \item \code{pi}: Matrix of size (\code{nrep}, \code{maxK}) with posterior samples of the distributional cluster weights.
  \item \code{omega}: Array of size (\code{maxL}, \code{maxK}, \code{nrep}) with observational cluster weights.
  \item \code{alpha}: Vector of length \code{nrep} with posterior samples of \eqn{\alpha}.
  \item \code{maxK}: Vector of length \code{nrep} with number of active distributional components.
}
}
\description{
\code{fit_fiSAN} fits the finite-infinite shared atoms nested (fiSAN) mixture model with Gaussian kernels and normal-inverse gamma priors on the unknown means and variances.
The function returns an object of class \code{SANmcmc} or \code{SANvi} depending on the chosen computational approach (MCMC or VI).
}
\details{
\strong{Data structure}

The finite-infinite common atoms mixture model is used to perform inference in nested settings, where the data are organized into \eqn{J} groups.
The data should be continuous observations \eqn{(Y_1,\dots,Y_J)}, where each \eqn{Y_j = (y_{1,j},\dots,y_{n_j,j})}
contains the \eqn{n_j} observations from group \eqn{j}, for \eqn{j=1,\dots,J}.
The function takes as input the data as a numeric vector \code{y} in this concatenated form. Hence, \code{y} should be a vector of length
\eqn{n_1+\dots+n_J}. The \code{group} parameter is a numeric vector of the same size as \code{y}, indicating the group membership for each
individual observation.
Notice that with this specification, the observations in the same group need not be contiguous as long as the correspondence between the variables
\code{y} and \code{group} is maintained.

\strong{Model}

The data are modeled using a Gaussian likelihood, where both the mean and the variance are observational-cluster-specific:
\deqn{y_{i,j}\mid M_{i,j} = l \sim N(\mu_l,\sigma^2_l)}
where \eqn{M_{i,j} \in \{1,\dots,L \}} is the observational cluster indicator of observation \eqn{i} in group \eqn{j}.
The prior on the model parameters is a normal-inverse gamma distribution \eqn{(\mu_l,\sigma^2_l)\sim NIG (m_0,\tau_0,\lambda_0,\gamma_0)},
i.e., \eqn{\mu_l\mid\sigma^2_l \sim N(m_0, \sigma^2_l / \tau_0)}, \eqn{1/\sigma^2_l \sim \text{Gamma}(\lambda_0, \gamma_0)} (shape, rate).

\strong{Clustering}

The model clusters both observations and groups.
The clustering of groups (distributional clustering) is provided by the allocation variables \eqn{S_j \in \{1,2,\dots\}}, with:
\deqn{Pr(S_j = k \mid \dots ) = \pi_k  \qquad \text{for } \: k = 1,2,\dots}
The distribution of the probabilities is \eqn{ \{\pi_k\}_{k=1}^{\infty} \sim GEM(\alpha)},
where GEM is the Griffiths-Engen-McCloskey distribution of parameter \eqn{\alpha},
which characterizes the stick-breaking construction of the DP (Sethuraman, 1994).

The clustering of observations (observational clustering) is provided by the allocation variables \eqn{M_{i,j} \in \{1,\dots,L\}}, with:
\deqn{ Pr(M_{i,j} = l \mid S_j = k, \dots ) = \omega_{l,k} \qquad \text{for } \: k = 1,2,\dots \, ; \: l = 1,\dots,L. }
The distribution of the probabilities is \eqn{(\omega_{1,k},\dots,\omega_{L,k})\sim \text{Dirichlet}_L(b,\dots,b)} for all \eqn{k = 1,2,\dots}.
Here, the dimension \eqn{L} is fixed.
}
\examples{
set.seed(123)
y <- c(rnorm(60), rnorm(40, 5))
g <- rep(1:2, rep(50, 2))
plot(density(y[g==1]), xlim = c(-5,10), main = "Group-specific density")
lines(density(y[g==2]), col = 2)

out_vi <- fit_fiSAN(y, group = g, est_method = "VI",
                    vi_param = list(n_runs = 1))
out_vi

out_mcmc <- fit_fiSAN(y = y, group = g, est_method = "MCMC")
out_mcmc
}

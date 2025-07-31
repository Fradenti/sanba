#' Compute and Plot the Posterior Similarity Matrix
#'
#' @description The function computes and plots the posterior similarity matrix (PSM), either for the whole dataset, or separately for each group.
#' The function takes as input an object from \code{fit_CAM}, \code{fit_fiSAN}, or \code{fit_fSAN},
#' used with the \code{est_method = "MCMC"} argument.
#'
#' @param object An object of class \code{SANmcmc}.
#' @param distributional Logical (default \code{FALSE}). If \code{FALSE}, the function computes the posterior similarity matrix (PSM) for the observational partition (i.e., between individual observations). If \code{TRUE}, it computes the PSM at the distributional level, that is, between groups.
#' @param group_specific Logical (default \code{FALSE}). If \code{FALSE}, the function considers the overall PSM.
#' If \code{TRUE}, the function considers the group-specific PSMs. This argument only affects the observational partition, i.e., when \code{distributional} is \code{FALSE}.
#' @param plot Logical (default \code{TRUE}). Whether to plot the PSM.
#' @param ncores A parameter to pass to the \code{salso::salso()} function. The number of CPU cores to use for parallel computing; a value of zero indicates the use of all cores of the system.
#'
#' @return The function \code{compute_psm} returns and plots the posterior similarity matrix.
#' When \code{distributional = FALSE}, if \code{group_specific = FALSE}, the output is a matrix of dimension \code{N x N};
#' if \code{group_specific = TRUE}, the output is a list on length \code{J} (the number of groups), where each entry contains a matrix of dimension \code{Nj x Nj}.
#' If \code{distributional = TRUE}, the output is a matrix of dimension \code{J x J}.
#'
#'
#'
#' @importFrom graphics image
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' y <- c(rnorm(100),rnorm(100,5))
#' g <- rep(1:2,rep(100,2))
#' plot(y,col=g)
#' # Fitting fiSAN via MCMC
#' est <- fit_fiSAN(y, g, est_method = "MCMC")
#' est
#' # Estimate PSM
#' psm_overall <- compute_psm(est)
#' # Estimate distributional PSM
#' psm_distrib <- compute_psm(est, distributional = TRUE)
#'
#' @export
#'
compute_psm <- function(object,
                        distributional = FALSE,
                        group_specific = FALSE, plot = TRUE, ncores = 0 ){
  if (!inherits(object, "SANmcmc")) {
    stop("compute_psm() is only defined for objects of class 'SANmcmc'.")
  }

  if(distributional){
    out <- salso::psm(object$sim$distr_cluster, nCores = ncores )
    if(plot) {
      n <- length(object$params$Nj)
      graphics::image(1:n, 1:n, out,
                      xlab = "Index", ylab = "Index",
                      main = "Posterior similarity matrix for all groups")
    }
  }else{

  if(!group_specific){
    out <- salso::psm(object$sim$obs_cluster, nCores = ncores )
    if(plot) {
      n <- length(object$params$y)
      graphics::image(1:n, 1:n, out,
            xlab = "Index", ylab = "Index",
            main = "Posterior similarity matrix for all observations")
    }

  }else{
    out <- list()
    J <- length(unique(object$params$group))
    for(j in 1:J){
      out[[j]] <- salso::psm(object$sim$obs_cluster[,object$params$group==j],
                            nCores = ncores)
    }

    if(plot){
      for(j in 1:J){
        nj <- object$params$Nj[j]
        graphics::image(1:nj, 1:nj, out[[j]],
              xlab = "Index", ylab = "Index", main = paste0("Posterior similarity matrix for group ", j))
        devAskNewPage(ask = TRUE)
      }
      devAskNewPage(ask = FALSE)
    }
  }
  }
  return(out)
}


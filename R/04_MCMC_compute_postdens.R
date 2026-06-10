#' Compute the Posterior Density for a Group
#'
#' @description
#' Computes posterior density estimates for a selected group from
#' a fitted SAN model estimated via MCMC. For each retained MCMC draw,
#' the function evaluates the corresponding posterior predictive density of
#' a new observation from the selected group.
#'
#' The returned object can be visualized with
#' \code{\link{plot.SANmcmc_postdens}}, which displays the collection of
#' posterior density curves together with pointwise posterior summaries.
#'
#' @param object An object of class \code{SANmcmc}, typically returned by
#'   \code{fit_CAM()}, \code{fit_fiSAN()}, or \code{fit_fSAN()} when
#'   \code{est_method = "MCMC"}.
#' @param group_ind Integer specifying the group for which posterior densities
#'   should be computed.
#' @param mcmc_considered Number of MCMC iterations (counting backward from the
#'   last saved iteration) used to compute posterior density samples.
#' @param lim Non-negative numeric value controlling the extension of the
#'   evaluation grid beyond the observed data range. The density is evaluated
#'   on
#'   \code{seq(min(y) - lim, max(y) + lim, length.out = length_yseq)}.
#' @param length_yseq Number of grid points used to evaluate the density.
#'
#' @return
#' An object of class \code{SANmcmc_postdens}, containing:
#' \describe{
#'   \item{densities}{
#'     A matrix of dimension \code{length_yseq x mcmc_considered}.
#'     Rows correspond to grid points and columns correspond to retained MCMC iterations.
#'   }
#'   \item{grid}{
#'     Numeric vector containing the evaluation grid.
#'   }
#'   \item{data}{
#'     Two-column matrix containing the observations and corresponding group
#'     labels for the selected group.
#'   }
#' }
#'
#' @seealso
#' \code{\link{fit_CAM}},
#' \code{\link{fit_fiSAN}},
#' \code{\link{fit_fSAN}},
#' \code{\link{plot.SANmcmc_postdens}}
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' y <- c(rnorm(100), rnorm(100, 5))
#' g <- rep(1:2, each = 100)
#'
#' # Fit fiSAN via MCMC
#' est <- fit_fiSAN(y, g, est_method = "MCMC")
#'
#' # Compute posterior density samples for group 1
#' postdens <- compute_postdens(
#'   est,
#'   group_ind = 2,
#'   mcmc_considered = 50
#' )
#'
#' # Plot posterior density summaries
#' plot(postdens)
#'
#' @export
compute_postdens <- function(object, group_ind = 1, mcmc_considered = 500,
                             lim = 1, length_yseq = 500){

  if (!inherits(object, "SANmcmc")) {
    stop("compute_postdens() is only defined for objects of class 'SANmcmc'.")
  }
  NSIM <- nrow(object$sim$mu)


  if (mcmc_considered > NSIM) {
    stop("'mcmc_considered' cannot exceed the number of saved MCMC iterations.")
  }

  if (!group_ind %in% unique(object$params$group)) {
    stop("Invalid group_ind value.")
  }


  # 1. Pre-allocate the grid
  x <- seq(min(object$params$y) - lim, max(object$params$y) + lim, length.out = length_yseq)

  # Determine dimensions
  maxL <- object$params$maxL

  # Pre-allocate the 3D array
  DENS <- array(NA, c(length(x),mcmc_considered))

  # 2. Vectorized computation over the iterations
  for (i in 1:mcmc_considered) {
    # Map back to the correct MCMC iteration index
    idx <- NSIM - i + 1

    # Extract components for the mixture
    mu    <- object$sim$mu[idx, 1:maxL]
    sigma <- sqrt(object$sim$sigma2[idx, 1:maxL])
    omega <- object$sim$omega[1:maxL, object$sim$distr_cluster[idx, group_ind], idx]

    # 3. Outer product magic:
    # outer() evaluates dnorm for every combination of x (rows) and mu/sigma (columns)
    # resulting in a length_yseq x maxL matrix.
    component_densities <- outer(x, 1:maxL, function(x_val, k) {
      dnorm(x_val, mean = mu[k], sd = sigma[k])
    })

    # Matrix multiplication (%*%) automatically multiplies by weights (omega)
    # and sums across the mixture components instantly.
    DENS[, i] <- component_densities %*% omega
  }

  whichg <- object$params$group[object$params$group==group_ind]
  whichy <- object$params$y[object$params$group==group_ind]

  D <- list(densities = DENS, grid=x,
            data = cbind(y = whichy, g = whichg))
  structure(D,
            class = c("SANmcmc_postdens", class(D)))

}

#' Plot Posterior Density Samples
#'
#' @description
#' Visualizes the output of \code{\link{compute_postdens}}. The plot displays
#' all posterior density curves obtained from the selected MCMC iterations,
#' together with pointwise posterior summaries and the observed data values.
#'
#' Specifically, the plot includes:
#' \itemize{
#'   \item all sampled density curves (gray lines),
#'   \item the posterior median density (black line),
#'   \item the posterior mean density (red line),
#'   \item pointwise credible interval bounds (blue lines),
#'   \item observed data values shown as rug-like points on the horizontal axis.
#' }
#'
#' @param x An object of class \code{SANmcmc_postdens} returned by
#'   \code{\link{compute_postdens}}.
#' @param alpha Significance level used to construct pointwise credible
#'   intervals. The default \code{alpha = 0.025} corresponds to a 95\%
#'   pointwise credible interval.
#' @param ... Additional graphical arguments passed to plotting methods.
#'
#' @return
#' Invisibly returns \code{NULL}. The function is called for its side effect of
#' producing a plot.
#'
#' @seealso
#' \code{\link{compute_postdens}}
#'
#' @importFrom graphics abline points lines matplot
#' @importFrom stats quantile
#'
#' @method plot SANmcmc_postdens
#' @export
plot.SANmcmc_postdens <- function(x, alpha = 0.025, ...){
  title <- paste0("Posterior distribution of group # ", unique(x$data[,2]))
  matplot(y = x$densities,x = x$grid, type="l", col="gray", lty=1,
          ylab = "Posterior density", xlab = "y", main = title, lwd = .5)
  graphics::abline(h=0,lty=2,lwd=.5)
  sum <- t(apply(x$densities, 1, function(xx) stats::quantile(xx, c(alpha, .5, 1-alpha))))
  graphics::lines(sum[,1]~x$grid,col=4)
  graphics::lines(sum[,2]~x$grid,col=1)
  graphics::lines(rowMeans(x$densities)~x$grid,col=2)
  graphics::lines(sum[,3]~x$grid,col=4)
  graphics::points(rep(0, nrow(x$data))~x$data[,1], cex = .5, pch=21, bg=4)
  }

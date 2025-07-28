#' Plot Variational Cluster Allocation Probabilities
#'
#' @description
#' Produces visualizations of the posterior cluster allocation probabilities from a SAN model fitted via variational inference.
#' The function supports plotting either the observation-level (OC) or distribution-level (DC) allocation probabilities, depending on the argument \code{distributional}.
#'
#' This function applies to objects returned by \code{\link{fit_CAM}}, \code{\link{fit_fiSAN}}, or \code{\link{fit_fSAN}} when used with \code{est_method = "VI"}.
#'
#' @param object An object of class \code{SANvi}, representing a model fitted via variational inference.
#' @param distributional Logical (default \code{FALSE}). If \code{FALSE}, plots the allocation probabilities of individual observations to observational clusters (OC). If \code{TRUE}, plots the allocation probabilities of groups to distributional clusters (DC).
#' @param ... Additional graphical parameters passed to the underlying \code{image()} function (or equivalent), for customizing the plot (e.g., \code{col}, \code{main}, \code{xlab}, \code{ylab}).
#'
#' @return The function plots the variational cluster allocation probabilities.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' y <- c(rnorm(60), rnorm(40, 5))
#' g <- rep(1:2, each = 50)
#'
#' # Fit fiSAN via VI
#' est <- fit_fiSAN(y, g, est_method = "VI")
#'
#' # Plot observational cluster probabilities
#' plot_vi_allocation_prob(est)
#'
#' # Plot distributional cluster probabilities
#' plot_vi_allocation_prob(est, distributional = TRUE)
#'
#' @export
plot_vi_allocation_prob <- function(object, distributional = FALSE, ...){
  if (!inherits(object, "SANvi")) {
    stop("plot_vi_allocation_prob() is only defined for objects of class 'SANvi'.")
  }
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))

  if(distributional){
    nrowRHO <- nrow(object$sim$RHO)
    ncolRHO <- ncol(object$sim$RHO)
    graphics::image(1:ncolRHO, 1:nrowRHO, t(object$sim$RHO), axes = F,
          xlab = "Cluster", ylab = "Group", main = "Distributional cluster allocation", ...)
    graphics::axis(1, c(0:ncolRHO)+.5, labels = NA)
    graphics::axis(1, seq(from = 0, to = ncolRHO, by=5), labels = seq(from=0, to=ncolRHO, by=5), col = NA, col.ticks = NA)
    graphics::axis(2, c(0:nrowRHO)+.5, labels = NA)
    graphics::axis(2, seq(from = 0, to = nrowRHO, by=1), labels = seq(from=0, to=nrowRHO, by=1), col = NA, col.ticks = NA)
    graphics::box()

  } else {

    graphics::par(mfrow= c(1,2))
    J <- object$params$J
    for(j in 1:J){
      nrowXI <- nrow(object$sim$XI[[j]])
      ncolXI <- ncol(object$sim$XI[[j]])
      graphics::image(1:ncolXI, 1:nrowXI, t(object$sim$XI[[j]]), axes = F,
                      xlab = "Cluster", ylab = "Observation", main = paste0("Group ", j), ...)
      graphics::axis(1, c(0:ncolXI)+.5, labels = NA)
      graphics::axis(1, seq(from = 0, to = ncolXI, by = 10), labels = seq(from = 0, to = ncolXI, by = 10), col = NA, col.ticks = NA)
      graphics::axis(2, c(0:nrowXI)+.5, labels = NA)
      graphics::axis(2, seq(from = 0, to = nrowXI, by = 10), labels = seq(from = 0, to = nrowXI, by = 10), col = NA, col.ticks = NA)
      graphics::box()
      devAskNewPage(ask = TRUE)
    }
    devAskNewPage(ask = FALSE)
  }
}

#' Estimate the Number of Observational and Distributional Clusters
#'
#' @description
#' Computes the estimated number of observational clusters (OC) and distributional clusters (DC) from a fitted SAN model object.
#'
#' For variational inference (\code{SANvi} objects), the function returns point estimates based on posterior mode assignments.
#' For MCMC-based inference (\code{SANmcmc} objects), it returns the mean, median, and variance of the number of clusters across posterior samples.
#'
#' @param object An object of class \code{SANvi} or \code{SANmcmc}, typically, the output of a call to \code{\link{fit_fiSAN}}, \code{\link{fit_fSAN}}, or \code{\link{fit_CAM}}.
#' @param ... ignored.
#'
#' @return A data frame reporting the estimated number of observational (OC) and distributional (DC) clusters.
#' \itemize{
#'   \item For \code{SANvi}: a single-row data frame with the point estimates.
#'   \item For \code{SANmcmc}: a three-row data frame with the mean, median, and variance across MCMC samples.
#' }
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' y <- c(rnorm(60), rnorm(40, 5))
#' g <- rep(1:2, each = 50)
#'
#' plot(density(y[g == 1]), xlim = c(-5, 10), main = "Group-specific density")
#' lines(density(y[g == 2]), col = 2)
#'
#' # Fit fiSAN via MCMC
#' est_mcmc <- fit_fiSAN(y, g, est_method = "MCMC")
#' number_clusters(est_mcmc)
#'
#' # Fit fiSAN via Variational Inference
#' est_vi <- fit_fiSAN(y, g, est_method = "VI")
#' number_clusters(est_vi)
#'
#' @rdname number_clusters
#' @export
number_clusters <- function(object, ...) {
  if (!inherits(object, "SANmcmc") && !inherits(object, "SANvi")) {
      stop("number_clusters() is only defined for objects of class 'SANmcmc' or 'SANvi'.")
        }
  UseMethod("number_clusters")
}

#' @export
number_clusters.SANmcmc <- function(object, ...){
  len_un <- apply(object$sim$obs_cluster, 1, function(x) length(unique(x)))

  len_un2 <- apply(object$sim$distr_cluster, 1, function(x) length(unique(x)))

  out <- data.frame( "Num1" = c(mean(len_un), median(len_un), var(len_un)), "Num2" = c(mean(len_un2), median(len_un2), var(len_un2)) )
  colnames(out) <- c("OC", "DC")
  rownames(out) <- c("Mean", "Median", "Variance")
  out
}

#' @export
number_clusters.SANvi <- function(object, ...){
  G <- length(unique(object$params$group))

  est_cl <- unlist( lapply(object$sim$XI,
                           function(X) apply(X, 1,
                                            which.max)))

  len_un <- length(unique(est_cl))

  est_cl2 <- apply(object$sim$RHO, 1, which.max)
  len_un2 <- length(unique(est_cl2))

  out <- data.frame( "OC" = len_un, "DC" = len_un2 )
  rownames(out) <- "Estimate"
  out
}



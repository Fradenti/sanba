#' Accessor functions for SAN model output
#'
#' @description
#' \code{get_model}, \code{get_time}, \code{get_params}, and \code{get_sim} are accessor functions that extract key components from objects of class \code{SANmcmc} or \code{SANvi},
#' depending on whether the model was estimated via MCMC or variational inference. In particular:
#'
#' \itemize{
#'   \item \code{get_model}: returns a character string indicating the model type;
#'   \item \code{get_time}: returns the runtime of the estimation procedure;
#'   \item \code{get_params}: returns the list of data and prior hyperparameters;
#'   \item \code{get_sim}: returns the fitted quantities (e.g., MCMC or VI estimates).
#' }
#'
#' @param object An object of class \code{SANmcmc} or \code{SANvi}.
#'
#' @return A component from the model object (see details).
#'
#' @examples
#' set.seed(123)
#' y <- c(rnorm(40, 0, 0.3), rnorm(20, 5, 0.3))
#' g <- c(rep(1:6, each = 10))
#' out <- fit_fSAN(y = y, group = g, est_method = "MCMC",
#'                 mcmc_param = list(nrep = 500, burn = 200))
#' get_model(out)
#' get_time(out)
#' hp <- get_params(out)
#' sims <- get_sim(out)
#'
#' @rdname get_accessors
#' @export
get_model <- function(object) {
  if (!inherits(object, "SANmcmc") && !inherits(object, "SANvi")) {
    stop("get_model() is only defined for objects of class 'SANmcmc' or 'SANvi'.")
  }
  UseMethod("get_model")
}

#' @rdname get_accessors
#' @export
get_model.SANvi <- function(object) {
  object$model
}

#' @rdname get_accessors
#' @export
get_model.SANmcmc <- function(object) {
  object$model
}

#' @rdname get_accessors
#' @export
get_time <- function(object) {
  if (!inherits(object, "SANmcmc") && !inherits(object, "SANvi")) {
    stop("get_time() is only defined for objects of class 'SANmcmc' or 'SANvi'.")
  }
  UseMethod("get_time")
}

#' @rdname get_accessors
#' @export
get_time.SANvi <- function(object) {
  object$time
}

#' @rdname get_accessors
#' @export
get_time.SANmcmc <- function(object) {
  object$time
}

#' @rdname get_accessors
#' @export
get_params <- function(object) {
  if (!inherits(object, "SANmcmc") && !inherits(object, "SANvi")) {
    stop("get_params() is only defined for objects of class 'SANmcmc' or 'SANvi'.")
  }
  UseMethod("get_params")
}

#' @rdname get_accessors
#' @export
get_params.SANvi <- function(object) {
  object$params
}

#' @rdname get_accessors
#' @export
get_params.SANmcmc <- function(object) {
  object$params
}

#' @rdname get_accessors
#' @export
get_sim <- function(object) {
  if (!inherits(object, "SANmcmc") && !inherits(object, "SANvi")) {
    stop("get_sim() is only defined for objects of class 'SANmcmc' or 'SANvi'.")
  }
  UseMethod("get_sim")
}

#' @rdname get_accessors
#' @export
get_sim.SANvi <- function(object) {
  object$sim
}

#' @rdname get_accessors
#' @export
get_sim.SANmcmc <- function(object) {
  object$sim
}

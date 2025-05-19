#' Print the variational inference output
#' @description Print method for objects of class \code{SANvi}.
#'
#' @param x Object of class \code{SANvi}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function prints a summary of the fitted model.
#'
#' @export
print.SANvi <- function(x, ... ){
  cat(paste("Variational inference results for", x$model ,"\n"))
  cat("----------------------------------------------\n")
  cat(paste("L:",x$params$L,"- K:",x$params$K,"\n"))
  cat(paste("Threshold:",x$params$epsilon,"\n"))
  cat(paste("ELBO value:", round(max(x$sim$Elbo_val),3),"\n"))
  cat(paste("Best run out of",x$params$n_runs,"\n"))
  cat(paste("Convergence reached in",length(x$sim$Elbo_val),"iterations\n"))
  cat(paste("Elapsed time:",round(as.numeric(x$time),3),units(x$time),"\n"))
  invisible(x)
}



#' Print the MCMC output
#' @description Print method for objects of class \code{SANmcmc}.
#'
#' @param x Object of class \code{SANmcmc}.
#' @param ... Ignored.
#'
#' @return The function prints a summary of the fitted model.
#'
#' @export
print.SANmcmc <- function(x, ...)
{
  cat("\n")
  cat(paste("MCMC result of", x$model, "model \n"))
  cat("-----------------------------------------------\n")
  cat(paste("Model estimated on", length(x$params$y), "total observations and",  length(unique(x$params$group)), "groups \n"))
  cat(paste("Size of the MCMC sample (after burn-in):", x$params$nrep - x$params$burn, "\n"))
  cat(paste("Total MCMC iterations performed:", x$params$nrep, "\n"))
  cat(paste("maxL:",x$params$maxL,"- maxK:",x$params$maxK,"\n"))
  cat(paste("Elapsed time:",round(as.numeric(x$time[[1]]),3), attr(x$time, "units"),"\n"))
}

#' Print the Variational Inference Output
#' @description Print method for objects of class \code{SANvi}.
#'
#' @param x Object of class \code{SANvi}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function prints a summary of the fitted model.
#'
#' @export
print.SANvi <- function(x, ... ){
  cat("\n")
  cat(paste("Variational inference results for", x$model ,"\n"))
  cat("----------------------------------------------\n")
  cat(paste("Model estimated on", length(x$params$y), "total observations and",  length(unique(x$params$group)), "groups \n"))
  if( length(x$params$Nj) <= 10 ){
    cat(paste("Groups sample sizes:", paste0(x$params$Nj, collapse = ", "), "\n\n"))
  }

  cat(paste("Threshold:",x$params$epsilon,"\n"))
  cat(paste("ELBO value:", round(max(x$sim$Elbo_val),3),"\n"))
  cat(paste("Best run out of",x$params$n_runs,"\n"))
  cat(paste("Convergence reached in",length(x$sim$Elbo_val),"iterations\n"))
  cat(paste("Elapsed time:",round(as.numeric(x$time),3),units(x$time),"\n"))
  cat("\n")
  invisible(x)
}



#' Print the MCMC Output
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
  cat(paste("MCMC results for", x$model, "\n"))
  cat("-----------------------------------------------\n")
  cat(paste("Model estimated on", length(x$params$y), "total observations and",  length(unique(x$params$group)), "groups \n"))
  if( length(x$params$Nj) <= 10 ){
  cat(paste("Groups sample sizes:", paste0(x$params$Nj, collapse = ", "), "\n\n"))
  }

  cat(paste("Size of the MCMC sample (after burn-in):", x$params$nrep - x$params$burn, "\n"))
  cat(paste("Total MCMC iterations performed:", x$params$nrep, "\n"))
  cat(paste("Elapsed time:",round(as.numeric(x$time[[1]]),3), attr(x$time, "units"),"\n"))
  cat("\n")
}

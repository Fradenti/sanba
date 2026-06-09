format_time <- function(x) {
  stopifnot(inherits(x, "difftime"))

  seconds <- as.numeric(x, units = "secs")

  h <- floor(seconds / 3600)
  m <- floor((seconds %% 3600) / 60)
  s <- seconds %% 60

  if (h > 0) {
    sprintf("%d h %d min %.1f s", h, m, s)
  } else if (m > 0) {
    sprintf("%d min %.1f s", m, s)
  } else {
    sprintf("%.3f s", s)
  }
}

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
  if(!is.null(x$all_difftimes)){
    cat(paste("Elapsed time (best run):",format_time(x$time),"\n"))
    cat(paste("Elapsed time (all runs):", format_time(Reduce("+",x$all_difftimes)),"\n"))
  }else{
    cat(paste("Elapsed time (single run):",format_time(x$time),"\n"))
  }
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
  cat(paste("Elapsed time:",format_time(x$time),"\n"))
  cat("\n")
}

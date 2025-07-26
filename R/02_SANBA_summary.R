#' Summarize the Variational Inference Output
#' @description Summary method for objects of class \code{SANvi}.
#'
#' @param object Object of class \code{SANvi}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function prints a summary of the fitted model.
#'
#' @export
summary.SANvi <- function(object, ...){
  G <- length(unique(object$params$group))

  est_cl <- unlist( lapply(object$sim$XI,
                           function(X) apply(X, 1,
                                             which.max)))

  len_un <- length(unique(est_cl))

  est_cl2 <- apply(object$sim$RHO, 1, which.max)
  len_un2 <- length(unique(est_cl2))

  out <- data.frame( "OC" = len_un, "DC" = len_un2 )
  rownames(out) <- "Estimate"


  cat("\n")
  cat(paste("Variational inference results for", object$model ,"\n"))
  cat("-----------------------------------------------\n")
  cat(paste("Model estimated on", length(object$params$y), "total observations and",  length(unique(object$params$group)), "groups \n"))

  if(length(object$params$Nj) <= 10){
    cat(paste("Groups sample sizes:", paste0(object$params$Nj, collapse = ", "), "\n\n"))
  }

  cat(paste("maxL:",object$params$L,"- maxK:",object$params$K,"\n"))
  cat(paste("Prior parameters (m0, tau0, lambda0, gamma0): (", paste0(c(object$params$m0, object$params$tau0,object$params$lambda0,object$params$gamma0), collapse = ", "), ") \n\n"))

  cat(paste("Threshold:",object$params$epsilon,"\n"))
  cat(paste("ELBO value:", round(max(object$sim$Elbo_val),3),"\n"))
  cat(paste("Best run out of",object$params$n_runs,"\n"))
  cat(paste("Convergence reached in",length(object$sim$Elbo_val),"iterations\n"))
  cat(paste("Elapsed time:",round(as.numeric(object$time),3),units(object$time),"\n\n"))

  cat("Number of observational and distributional clusters:\n")
  print(out)
  invisible(object)
}



#' Summarize the MCMC Output
#' @description Summary method for objects of class \code{SANmcmc}.
#'
#' @param object of class \code{SANmcmc}.
#' @param ... Ignored.
#'
#' @return The function prints a summary of the fitted model.
#'
#' @export
summary.SANmcmc <- function(object, ...)
{
  len_un <-  apply(object$sim$obs_cluster, 1, function(x) length(unique(x)))
  len_un2 <- apply(object$sim$distr_cluster, 1, function(x) length(unique(x)))

  out <- data.frame( "Num1" = c(mean(len_un), median(len_un), var(len_un)),
                    "Num2" = c(mean(len_un2), median(len_un2), var(len_un2)) )
  colnames(out) <- c("OC", "DC")
  rownames(out) <- c("Mean", "Median", "Variance")


  cat("\n")
  cat(paste("MCMC results for", object$model, "\n"))
  cat("-----------------------------------------------\n")
  cat(paste("Model estimated on", length(object$params$y), "total observations and",  length(unique(object$params$group)), "groups \n"))
  if( length(object$params$Nj) <= 10 ){
    cat(paste("Groups sample sizes:", paste0(object$params$Nj, collapse = ", "), "\n\n"))
  }
  cat(paste("maxL:",object$params$maxL,"- maxK:",object$params$maxK,"\n"))
  cat(paste("Prior parameters (m0, tau0, lambda0, gamma0): (", paste0(c(object$params$m0, object$params$tau0,object$params$lambda0,object$params$gamma0), collapse = ", "), ") \n\n"))

  cat(paste("Size of the MCMC sample (after burn-in):", object$params$nrep - object$params$burn, "\n"))
  cat(paste("Total MCMC iterations performed:", object$params$nrep, "\n"))
  cat(paste("Elapsed time:",round(as.numeric(object$time[[1]]),3), attr(object$time, "units"),"\n\n"))

  cat("Number of observational and distributional clusters:\n")
  print(out)
  invisible(object)
}


#' Compute posterior mean of stick-breaking weights.
#'
#' @param alk,blk the posterior variational parameters for the SB betas.
#' @noRd
#'
post_sb_weight <- function(alk, blk){
  
  K  <- ncol(alk)
  L  <- nrow(alk)
  # are we sure about this?
  logomega_post <- matrix(NA,L,K)
  
  for(k in 1:K){
    
    p2 <- blk[,k]/(alk[,k]+blk[,k])
    p2 <- c(1,p2[-L])
    
    logomega_post[,k] <- log(alk[,k]/(alk[,k]+blk[,k])) + cumsum(log(p2))
    
  }
  
  return(exp(logomega_post))
  
}



#' Relabel clusters 
#' @keywords internal
#' @param ix vector of cluster labels
.relabel <- function(ix) 
{
  if(min(ix)==0) ix <- ix+1
  while( max(ix) != length(unique(ix)) )
  {
    missing_label1 <- (1:max(ix))[sapply(1:max(ix), function(x) !(x %in% ix) )][1]
    ix[ix>missing_label1] <- ix[ix>missing_label1]-1
  }
  return(ix)
}


#' Extract best
#' @param object 
#' @keywords internal
extract_best <- function(object) {

  summa <- do.call(rbind,
                   lapply(object,  function(q)
                     c(max = max(q$sim$Elbo_val))))
  ind <- which.max(summa)

  return(object[[ind]])
}


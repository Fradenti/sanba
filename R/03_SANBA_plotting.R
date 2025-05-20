#' Visual check of convergence of the MCMC output
#' @description Plot method for objects of class \code{SANmcmc}.
#' Check the convergence of the MCMC through visual inspection of the chains.
#'
#' @param x Object of class \code{SANmcmc} (usually, the result of a call to \code{fit_CAM}, \code{fit_fiSAN}, or \code{fit_fSAN}, used with the \code{est_method = "MCMC"} argument).
#' @param param String with the names of the parameters to check. It can be one of \code{"mu"}, \code{"sigma2"}, \code{"pi"},
#' \code{"num_clust"}, \code{"alpha"}, \code{"beta"}.
#' @param show_density  Logical (default \code{TRUE}). Should a kernel estimate of the density be plotted?
#' @param add_burnin Integer (default = 0). Additional number of observations to discard in the burn-in.
#' @param show_convergence Logical (default \code{TRUE}). Should a superimposed red line of the cumulative mean be plotted?
#' @param trunc_plot Integer (default = 10). For multidimensional parameters, the maximum number of components to be plotted.
#' @param ... Ignored.
#' @note The function is not available for the observational weights \eqn{\omega}.
#'
#' @return The function displays the traceplots and posterior density estimates of the parameters sampled in the MCMC algorithm.
#'
#' @examples
#' set.seed(123)
#' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' g <- c(rep(1,30), rep(2, 30))
#' out <- fit_fiSAN(y = y, group = g, "MCMC", mcmc_param = list(nrep = 500, burn = 200))
#' plot(out, param = "mu", trunc_plot = 2)
#' plot(out, param = "sigma2", trunc_plot = 2)
#' plot(out, param = "alpha", trunc_plot = 1)
#' plot(out, param = "alpha", add_burnin = 100)
#' plot(out, param = "pi", trunc_plot = 4, show_density = FALSE)
#'
#' out <- fit_CAM(y = y, group = g, "MCMC",
#' mcmc_param = list(nrep = 500, burn = 200, seed= 1234))
#' plot(out, param = "mu", trunc_plot = 2)
#' plot(out, param = "sigma2", trunc_plot = 2)
#' plot(out, param = "alpha")
#' plot(out, param = "pi", trunc_plot = 2)
#' plot(out, param = "pi", trunc_plot = 5)
#' plot(out, param = "num_clust", trunc_plot = 5)
#' plot(out, param = "beta", trunc_plot = 2)
#'
#' out <- fit_fSAN(y = y, group = g, "MCMC", mcmc_param = list(nrep = 500, burn = 200))
#' plot(out, param = "mu", trunc_plot = 2)
#' plot(out, param = "sigma2", trunc_plot = 2)
#' plot(out, param = "pi", trunc_plot = 4,
#'      show_convergence = FALSE, show_density = FALSE)
#'
#' @importFrom graphics par
#' @importFrom grDevices devAskNewPage
#' @export
plot.SANmcmc <- function(x, param = c("mu",
                                       "sigma2",
                                       "pi",
                                       "num_clust",
                                       "alpha",
                                       "beta"),
                         show_density = TRUE,
                         add_burnin = 0,
                         show_convergence = TRUE,
                         trunc_plot = 2,
                         ...)
{


  param <- match.arg(param)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if( (x$model == "fiSAN") & param == "beta"){
    stop("beta is not available for fiSAN models")
  }
  if( (x$model == "fSAN") ){
    if(param == "alpha" | param == "beta"){
    stop("alpha or beta are not available for fSAN models")
  }
  }

  if(show_density){
    .traces_and_density(x$sim, param,
            add_burnin,
            show_convergence,
            trunc_plot)
  } else {
    .traces(x$sim, param,
            add_burnin,
            show_convergence,
            trunc_plot)
  }
  devAskNewPage(ask = F)
}



#' @importFrom graphics par
#' @importFrom grDevices devAskNewPage
#' @keywords internal
.traces <- function(sim, param,
                    add_burnin,
                    show_convergence,
                    trunc_plot)
{


  count <- 0
  tmp_names <- c()


  if(param == "num_clust"){
    udc <- apply(sim$distr_cluster,1,function(x) length(unique(x)))
    uoc <- apply(sim$obs_cluster,1,function(x) length(unique(x)))
    tmp <- cbind("# DCs" = udc, "# OCs" = uoc)
    cl_names <- colnames(tmp)
  }else{
    stringg <- paste0("sim$", param)
    tmp <- eval(parse(text=stringg))
  }
  if(dim(tmp)[2] > trunc_plot) { tmp_names <- c(tmp_names, param) }
  for(j in 1:min(dim(tmp)[2],  trunc_plot)){
    count <- count+1
  }

  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))
  graphics::par(mfrow= c(min(3,count),2), mar = c(4, 4, 2, 2) + 0.1)

  devAskNewPage(ask <- F)
  for(j in 1:min(dim(tmp)[2],  trunc_plot)){

    if(param == "num_clust"){
      namem <- cl_names[j]
    }else if(dim(tmp)[2] == 1) { namem <- paste0(param) } else { namem <- paste0(param,"_",j) }

    if(add_burnin>0){
      plot(tmp[-(1:add_burnin),j], type="l", main = namem , xlab = "Iteration", ylab = "Value")
      if(show_convergence) {
        lines(1:length(tmp[-(1:add_burnin),j]), cumsum((tmp[-(1:add_burnin),j]))/(1:length(tmp[-(1:add_burnin),j])), col=2)
      }
    }else{
      plot(tmp[,j], type="l", main = namem , xlab = "Iteration", ylab = "Value")
      if(show_convergence) {
        lines(1:length(tmp[,j]), cumsum((tmp[,j]))/(1:length(tmp[,j])), col=2)
      }

    }

    devAskNewPage(ask = T)
  }
  if(dim(tmp)[2]>trunc_plot){
    print(paste0("Output truncated at ", trunc_plot, " for ", param, "."))
  }else{
    print(paste0("Output for ", param, "."))
  }
}




#' @importFrom graphics par
#' @importFrom grDevices devAskNewPage
#' @importFrom stats density
#' @keywords internal
.traces_and_density <- function(sim, param,
                                add_burnin,
                                show_convergence,
                                trunc_plot)
{


  if("omega" %in% param) {
    warning("Traceplots for omega are not available")
    param <- param[param!="omega"]
  }

  count <- 0
  tmp_names <- c()


  if(param == "num_clust"){
    udc <- apply(sim$distr_cluster,1,function(x) length(unique(x)))
    uoc <- apply(sim$obs_cluster,1,function(x) length(unique(x)))
    tmp <- cbind("# DCs" = udc, "# OCs" = uoc)
    cl_names <- colnames(tmp)
  }else{
    stringg <- paste0("sim$", param)
    tmp <- eval(parse(text=stringg))
  }
  if(dim(tmp)[2] > trunc_plot) { tmp_names <- c(tmp_names, param) }
    for(j in 1:min(dim(tmp)[2],  trunc_plot)){
      count <- count+1
    }

  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))
  graphics::par(mfrow= c(min(3,count),2), mar = c(4, 4, 2, 2) + 0.1)

  devAskNewPage(ask = F)

    for(j in 1:min(dim(tmp)[2],  trunc_plot)){

      if(param == "num_clust"){
        namem <- cl_names[j]
      }else if(dim(tmp)[2] == 1) { namem <- paste0(param) } else { namem <- paste0(param,"_",j) }

      if(add_burnin>0){
      plot(tmp[-(1:add_burnin),j], type="l", main = namem , xlab = "Iteration", ylab = "Value")
      if(show_convergence) {
        lines(1:length(tmp[-(1:add_burnin),j]), cumsum((tmp[-(1:add_burnin),j]))/(1:length(tmp[-(1:add_burnin),j])), col=2)
        if(param == "num_clust"){
          plot(table(tmp[,j]), main = namem)
        }else{
          plot(density(tmp[,j]), main = namem)
        }
      }
      }else{
        plot(tmp[,j], type="l", main = namem , xlab = "Iteration", ylab = "Value")
        if(show_convergence) {
          lines(1:length(tmp[,j]), cumsum((tmp[,j]))/(1:length(tmp[,j])), col=2)

          if(param == "num_clust"){
            plot(table(tmp[,j]), main = namem)
          }else{
              plot(density(tmp[,j]), main = namem)
          }
        }

      }

      devAskNewPage(ask = T)
    }

  if(dim(tmp)[2]>trunc_plot){
    print(paste0("Output truncated at ", trunc_plot, " for ", param, "."))
  }else{
    print(paste0("Output for ", param, "."))
  }
}

#' Visual check of convergence of the VI output
#'
#' @description Plot method for objects of class \code{SANvi}.
#' The function displays two graphs. The left plot shows the progression of all the ELBO values as a function of the iterations.
#' The right plots shows the ELBO increments between successive iterations of the best run on a log scale (note: increments should always be positive).
#'
#' @param x Object of class \code{SANvi} (usually, the result of a call to \code{fit_CAM}, \code{fit_fiSAN}, or \code{fit_fSAN}, used with the \code{est_method = "VI"} argument).
#'
#' @param ... Ignored.
#'
#' @return The function plots the path followed by the ELBO and its subsequent differences.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- c(rnorm(200,0,0.3), rnorm(100,5,0.3))
#' g <- c(rep(1,150), rep(2, 150))
#' out <- fit_fSAN(y = y, group = g, "VI", vi_param = list(n_runs = 2))
#' plot(out)
plot.SANvi <- function(x, ...){

  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))
  graphics::par(mfrow= c(1,2))
  if(is.null(x$all_elbos)){
    plot(x$sim$Elbo_val,
       xlab = "Iterations - log scale", ylab = "ELBO",
       main = paste(x$model, "- ELBO"),type="b",cex=.5,
       log="x")
  }else{
      lli = lapply(x$all_elbos, length)
      lmi = lapply(x$all_elbos, min)
      lma = lapply(x$all_elbos, max)
      lims = c(min(unlist(lmi)),max(unlist(lma)))

      plot(x$all_elbos[[1]],
           xlab = "Iterations - log scale", ylab = "ELBO",
           main = paste(x$model, "- Results over",length(lli),"runs\nELBO trajectories"), ylim = lims,
           xlim = c(1,max(unlist(lli))),
           type="b",cex=.5,
           log="x", col = "gray")
    for(i in 1:length(x$all_elbos)){
      points(x$all_elbos[[i]],
           xlab = "Iterations - log scale", ylab = "ELBO",
           type="b",cex=.5, col = "gray")
    }
      points(x$sim$Elbo_val,
             xlab = "Iterations - log scale", ylab = "ELBO",
             type="b",cex=.5, col = "steelblue3")
  }
  plot(diff(x$sim$Elbo_val),
       xlab = "Iterations - log scale", ylab = "diff(ELBO)",
       main = paste(x$model, "- Best run\nELBO differences"),type="b",cex=.5,
       log="x")

  }






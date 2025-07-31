#' Estimate the Atoms and Weights of the Discrete Mixing Distributions
#'
#' The function computes the posterior means of the atoms and weights characterizing the discrete mixing distributions.
#' The function takes as input an object from \code{fit_CAM}, \code{fit_fiSAN},
#' or \code{fit_fSAN}, used with the \code{est_method = "VI"} argument, and returns an object of class \code{SANvi_G}.
#'
#' @param object An object of class \code{SANvi}.
#'
#' @return The function \code{estimate_G} returns an object of class \code{SANvi_G}, which is a matrix comprising the posterior means,
#' variances, and weights of each estimated DC (one mixture component for each row).
#'
#' @export
#'
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' y <- c(rnorm(100),rnorm(100,5))
#' g <- rep(1:2,rep(100,2))
#' plot(y,col=g)
#' # Fitting fiSAN via variational inference
#' est <- fit_fiSAN(y,g,vi_param= list(n_runs = 10))
#' est
#' summary(est)
#' # Estimate posterior atoms and weights
#' G <- estimate_G(est)
#' summary(G)
estimate_G <- function(object){

  if (!inherits(object, "SANvi")) {
    stop("estimate_G() can only be applied to objects of class 'SANvi'.")
  }

  if(object$model == "CAM"){
    exp_weights <- post_sb_weight(object$sim$a_bar_lk, object$sim$b_bar_lk)

  }else   if(object$model == "fSAN" |
             object$model == "fiSAN"){
    exp_weights <- apply(object$sim$b_dirichlet_lk,2,function(x) x/sum(x))
  }else{
    stop("Provide valid variational inference object")
  }

  exp_sigma2 <- object$sim$theta_l[,4]/(object$sim$theta_l[,3]-1)
  exp_mu     <- object$sim$theta_l[,1]
  ind_row  <-  unique(unlist(lapply(object$sim$XI,
                                    function(x) apply(x, 1, which.max))))
  ind_col <- unique(apply(object$sim$RHO,1,which.max))


  W <- exp_weights[ind_row,ind_col]

  D <- data.frame(post_mean = exp_mu[ind_row],
                  post_var  = exp_sigma2[ind_row],
                  post_weigths_DC = W)


  if(is.null(nrow(W))){
    colnames(D)[-c(1:2)] <- paste0("post_weigths_DC",1:(ncol(D)-2))
  }else{
    wmeans <- apply(W,2,function(x) sum(x*D$post_mean))
    D <- D[order(D$post_mean),c(1:2,order(wmeans)+2)]
    colnames(D)[-c(1:2)] <- paste0("post_weigths_DC",1:(ncol(D)-2))
  }

  structure(D,
            class = c("SANvi_G", class(D)))

}

#' @name estimate_G
#'
#' @param x an object of class \code{summary_SANvi} (usually, the result of a call to \code{fit_CAM}, \code{fit_fiSAN}, or \code{fit_fSAN}, used with the \code{est_method = "VI"} argument).
#' @param DC_num an integer or a vector of integers indicating which distributional clusters to plot.
#' @param lim optional value for the \code{plot} method to adjust the limits of the x-axis (the default is 2). The atoms are plotted on a range
#' given by \code{min(posterior means)-lim, max(posterior means)+lim}.
#'
#' @importFrom graphics abline points lines
#' @importFrom stats dnorm
#'
#' @export
#'
plot.SANvi_G <- function(x,
                                  DC_num = NULL,
                                  lim = 2, ...) {

  if (is.null(DC_num)) {
    DC_num <- 1:(ncol(x)-2)
  }
  if(max(DC_num) > (ncol(x)-2) ){
    stop(paste0("There are less estimated DCs than requested.\n",
                "Please provide a number for DC_num between 1 and ", (ncol(x)-2),"." ))
  }

  atoms   <-  x[, 1:2]
  seqq    <- seq(min(atoms[, 1]) - lim,
                 max(atoms[, 1]) + lim,
                 length.out = 250)

  width <- sqrt(atoms[, 2])

  weights <-  x[, DC_num[1] + 2]
  dens <- sapply(seqq, function(g) {
    sum(weights * stats::dnorm(g, atoms[, 1], sqrt(atoms[, 2])))
  })
  norm_dens <- dens / max(dens)

  if(length(DC_num) == ncol(x)-2){
    main_title <- paste0("All distributional clusters")
  }else if(length(DC_num) > 8){
    main_title <- paste0("Distributional clusters")
  }else{
    main_title <- paste0("Distributional clusters # ", paste(DC_num, collapse = ", "))
  }

  plot(
    weights ~ atoms[, 1],
    type = "h",
    lwd = 2,
    lend = 1,
    xlim = c(min(seqq), max(seqq)),
    ylim = c(0, 1),
    xlab = "y",
    ylab = "Posterior weights and normalized density",
    main = main_title
  )
  graphics::points(weights ~ atoms[, 1], cex = width)
  graphics::lines(norm_dens ~ seqq)

  if (length(DC_num) > 1) {
    for (j in 2:length(DC_num)) {
      weights <-  x[, DC_num[j] + 2]
      dens <- sapply(seqq, function(g) {
        sum(weights * stats::dnorm(g, atoms[, 1], sqrt(atoms[, 2])))
      })
      norm_dens <- dens / max(dens)

      graphics::points(
        weights ~ atoms[, 1],
        type = "h",
        lwd = 2,
        lend = 1,
        xlim = c(min(seqq), max(seqq)),
        ylim = c(0, 1),
        col = j
      )
      graphics::points(weights ~ atoms[, 1], cex = width, col = j)
      graphics::lines(norm_dens ~ seqq, col = j)
    }

  }

}

#' @name estimate_G
#'
#' @param object an object of class \code{SANvi_G} (usually, the result of a call to \code{estimate_G}).
#' @param thr argument for the \code{print} method. It should be a small positive number,
#' representing a threshold. If the posterior weight of a specific shared atom is below the threshold, the
#' atom is not reported.
#' @param ... ignored.
#'
#' @export
#'
summary.SANvi_G <- function(object, thr = 1e-2, ...){
  cat(paste("Estimated random measures via VI\n"))
  cat("----------------------------------\n")

  rows <- list()
  for(i in 1:(ncol(object)-2)){
    rows[[i]] <- which(object[,2+i] > thr)
  }
  n_g <- length(rows)

  cat(paste("Atoms with posterior weight >", thr, "\n"))
  # cat("----------------------------------\n")
  cat(paste("Number of detected DCs:", n_g, "\n"))
  cat("----------------------------------\n")

  for(j in 1:n_g){

    cat(paste("\nDistributional cluster #",j,"\n"))

    if(length(rows[[j]])==0){
      cat(paste0("No atom has weight above the selected threshold of ",thr,"\n"))
      next()
    }
    Dsubj <- round(object[rows[[j]],c(1,2,j+2)],3)

    if(is.null(nrow(Dsubj))){
      Dsubj <- matrix(Dsubj,1,3)
      rownames(Dsubj) <- "1"
      colnames(Dsubj) <- c("post_mean", "post_var", "post_weight")
    }else{
      colnames(Dsubj)[3] <- "post_weight"
    }
    print(data.frame(Dsubj))

  }
}


#' @name estimate_G
#'
#' @param x an object of class \code{SANvi_G} (usually, the result of a call to \code{estimate_G}).
#' @param thr argument for the \code{print} method. It should be a small positive number,
#' representing a threshold. If the posterior weight of a specific shared atom is below the threshold, the
#' atom is not reported.
#' @param ... ignored.
#'
#' @export
#'
print.SANvi_G <- function(x, thr = 1e-2, ...){
  G2 <- as.data.frame(x)
  atoms <- G2[,1:2]
  print(
    cbind(round(atoms,3),
    as.data.frame(
      apply(round(as.matrix(G2[,-c(1:2)]),3), 2 ,
            function(z) ifelse(z<thr, ".",z)))))
}

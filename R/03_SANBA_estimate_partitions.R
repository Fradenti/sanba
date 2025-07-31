#' Estimate the Observational and Distributional Partition
#'
#' @description Given the output of a \code{sanba} model-fitting function, this method estimates both the observational and distributional partitions.
#' For MCMC objects, it computes a point estimate using \code{\link[salso:salso]{salso::salso()}};
#' for Variational Inference (VI) objects, the cluster allocation is determined by the label with the highest estimated variational probability.
#'
#' @param object Object of class \code{SANmcmc} (usually, the result of a call to \code{\link{fit_fiSAN}},
#' \code{\link{fit_fSAN}}, or \code{\link{fit_CAM}} with \code{est_method = "MCMC"}) or \code{SANvi}
#' (the result of a call to \code{\link{fit_fiSAN}}, \code{\link{fit_fSAN}}, or \code{\link{fit_CAM}} with \code{est_method = "VI"}).
#' @param add_burnin Integer (default = 0). Number of observations to discard as additional burn-in (only for \code{SANmcmc} objects).
#' @param ncores A parameter to pass to the \code{salso::salso()} function (only for \code{SANmcmc} objects). The number of CPU cores to use for parallel computing; a value of zero indicates the use of all cores on the system.
#' @param ordered Logical, if \code{TRUE} (default), the function sorts the distributional cluster labels reflecting the
#' increasing values of medians of the data assigned to each DC. If \code{FALSE}, no ordering is applied.

#' @param ... further arguments passed to or from other methods.
#'
#' @return A list of class \code{partition_vi} or \code{partition_mcmc} containing
#' \itemize{
#'   \item \code{obs_level}: a data frame containing the data values, their group indexes, and the observational and distributional clustering assignments for each observation.
#'   \item \code{dis_level}: a vector with the distributional clustering assignment for each unit.
#' }
#'
#'
#' @rdname estimate_partition
#'
#' @export
estimate_partition <- function(object, ...) {
  if (!inherits(object, "SANmcmc") && !inherits(object, "SANvi")) {
    stop("get_model() is only defined for objects of class 'SANmcmc' or 'SANvi'.")
  }
  UseMethod("estimate_partition")
}

#' @rdname estimate_partition
#' @examples
#' set.seed(123)
#' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' g <- c(rep(1:6, each = 10))
#' out <- fit_fSAN(y = y, group = g, "VI", vi_param = list(n_runs = 10))
#' plot(out)
#' clust <- estimate_partition(out)
#' summary(clust)
#' plot(clust, lwd = 2, alt_palette = TRUE)
#' plot(clust, type = "scatter", alt_palette = FALSE, cex = 2)
#'
#' @export
estimate_partition.SANvi <- function(object, ordered = TRUE, ...) {

  if(!inherits(object, "SANvi")){
    warning("The passed object should be of class 'SANvi'")
  }

  estimated_oc <- as.numeric(factor(unlist(lapply(object$sim$XI,
                                               function(x) apply(x, 1, which.max)))))
  estimated_dc <- as.numeric(factor(apply(object$sim$RHO,1,which.max)))

  if(ordered){
    estimated_oc <- unname(rank(tapply(object$params$y,
                                       INDEX = estimated_oc,
                                       median))[estimated_oc])

    estimated_dc <- unname(rank(tapply(object$params$y,
                                       INDEX = rep(estimated_dc,object$params$Nj),
                                       median))[estimated_dc])
  }
  D <- data.frame(Y  = object$params$y,
                  G  = object$params$group,
                  OC = estimated_oc,
                  DC = rep(estimated_dc,object$params$Nj))

  L <- (list(obs_level = D,
             dis_level = estimated_dc))

  structure(L,class = c("partition_vi",class(L)))
  }


#' @rdname estimate_partition
#' @seealso \code{\link[salso:salso]{salso::salso()}}
#'
#' @importFrom salso salso
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
#' g <- c(rep(1:6, each = 10))
#' out <- fit_fSAN(y = y, group = g, "MCMC", mcmc_param=list(nrep=500,burn=200))
#' plot(out)
#' clust <- estimate_partition(out)
#' summary(clust)
#' plot(clust, lwd = 2)
#' plot(clust,  type = "boxplot", alt_palette = TRUE)
#' plot(clust,  type = "scatter", alt_palette = TRUE, cex = 2, pch = 4)
#'
estimate_partition.SANmcmc <- function(object, ordered = TRUE, add_burnin = 0, ncores = 0, ...) {


  if(!inherits(object, "SANmcmc")){
    warning("The passed object should be of class 'SANmcmc'")
  }
  if(add_burnin>0) {
    OC <- object$sim$obs_cluster[-1:add_burnin,]
    DC <- object$sim$distr_cluster[-1:add_burnin,]
  }else{
    OC <- object$sim$obs_cluster
    DC <- object$sim$distr_cluster
  }

  estimated_oc <- suppressWarnings(salso::salso(OC, nCores = ncores))
  estimated_dc <- suppressWarnings(salso::salso(DC, nCores = ncores))

  if(ordered){
    estimated_oc <- unname(rank(tapply(object$params$y,
                                    INDEX = estimated_oc,
                                    median))[estimated_oc])

    estimated_dc <- unname(rank(tapply(object$params$y,
                                     INDEX = rep(estimated_dc,object$params$Nj),
                                     median))[estimated_dc])
  }
  D <- data.frame(Y  = object$params$y,
                  G  = object$params$group,
                  OC = estimated_oc,
                  DC = rep(estimated_dc,object$params$Nj))

  L <- (list(obs_level = D,
             dis_level = estimated_dc))

  structure(L, class = c("partition_mcmc", class(L)))
}

#' @name estimate_partition
#'
#' @export
#'
summary.partition_mcmc <- function(object, ...){
  cat("Summary of the posterior observ. and distrib. partitions estimated via MCMC\n")
  cat("----------------------------------\n")
  cat(paste("Number of estimated OCs:", length(unique(object$obs_level$OC)),"\n"))
  cat(paste("Number of estimated DCs:",length(unique(object$dis_level)),"\n"))
  cat("----------------------------------\n")
  invisible(object)
}


#' @name estimate_partition
#'
#' @export
#'
summary.partition_vi <- function(object, ...){
  cat("Summary of the posterior observ. and distrib. partitions estimated via VI\n")
  cat("----------------------------------\n")
  cat(paste("Number of estimated OCs:", length(unique(object$obs_level$OC)),"\n"))
  cat(paste("Number of estimated DCs:",length(unique(object$dis_level)),"\n"))
  cat("----------------------------------\n")
  invisible(object)
}

#' @name estimate_partition
#'
#' @export
#'
print.partition_mcmc <- function(x, ...){
  cat("Estimated posterior partitions via MCMC\n")
  invisible(x)
}


#' @name estimate_partition
#'
#' @export
#'
print.partition_vi <- function(x, ...){
  cat("Estimated posterior partitions via VI\n")
  invisible(x)
}


#' @name estimate_partition
#'
#' @param x The result of a call to \code{\link{estimate_partition}}.
#' @param DC_num An integer or a vector of integers indicating which distributional clusters to plot.
#' @param type What type of plot should be drawn. Available types are "boxplot", "ecdf", and "scatter".
#' @param alt_palette Logical, the color palette to be used. Default is \code{R} base colors (\code{alt_palette = FALSE}).
#' @param ... Additional graphical parameters to be passed to the \code{plot} function.
#'
#' @importFrom graphics abline lines points boxplot
#' @importFrom stats median
#' @importFrom grDevices colorRampPalette
#' @importFrom scales alpha
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
plot.partition_mcmc <- function(x,
                         DC_num = NULL,
                         type = c("ecdf", "boxplot", "scatter"),
                         alt_palette = FALSE,
                         ...) {

  type   <- match.arg(type)
  if (is.null(DC_num)) {
    DC_num <- 1:max(x$dis_level)
  }
  if(max(DC_num) > max(x$dis_level)){
    stop(paste0("There are less estimated DCs than requested.\n",
                "Please provide a number for DC_num between 1 and ", max(x$dis_level) ))
  }

  max_OC <- max(x$obs_level$OC)
  max_DC <- max(x$dis_level)
  if(alt_palette){
    colpal_DC <- colorRampPalette(brewer.pal(8, "Dark2"))(max_DC)
    colpal_OC <- colorRampPalette(brewer.pal(8, "Dark2"))(max_OC)
  }else{
    colpal_DC <- 1:max_DC
    colpal_OC <- 1:max_OC
  }


  dix <- rank(tapply(x$obs_level$Y, x$obs_level$DC, stats::median))

  ind_ord_dis <- sapply(x$dis_level, function(g) dix[g])
  ind_ord_obs <- dix[x$obs_level$DC]


  inds_row <- which(ind_ord_obs %in% DC_num)
  inds_col <- which(ind_ord_dis %in% DC_num)

  suby <- (x$obs_level$Y[inds_row])
  subg <- (x$obs_level$G[inds_row])
  subDC <- x$obs_level$DC[inds_row]
  subOC <- x$obs_level$OC[inds_row]

  if(length(DC_num) == max(x$dis_level)){
    main_title <- paste0("All distributional clusters")
  }else if( length(DC_num) > 8 ){
    main_title <- NULL
  }else{
    main_title <- paste0("Distributional clusters # ", paste(DC_num, collapse = ", "))
  }



  if (type == "ecdf") {
    X <- unique(x$obs_level$G[inds_row])
    Nj <- table(x$obs_level$G)
    ysteps <- (0:(Nj[inds_col][1] - 1)) /
      (Nj[inds_col][1] - 1)
    xsteps <- sort(suby[subg == X[1]])


    plot(
      ysteps ~ xsteps,
      type = "l",
      xlim = c(min(suby), max(suby)),
      ylim = c(0, 1),
      col = scales::alpha(colpal_DC[ind_ord_dis[inds_col][1]], .5),
      xlab = "y",
      ylab = "eCDF",
      main = paste0("eCDFs colored by DC\n", main_title), ...
    )
    graphics::points(
      (ysteps ~ xsteps),
      cex = .1,
      col = scales::alpha(colpal_DC[ind_ord_dis[inds_col][1]], .5)
    )
    graphics::abline(h = c(0, 1),
                     col = "gray",
                     lty = 3)


    for (j in 2:length(X)) {
      ysteps <- (0:(Nj[X[j]] - 1)) / (Nj[X[j]] - 1)
      xsteps <- sort(suby[subg == X[j]])

      graphics::lines((ysteps ~ xsteps),
                      col = scales::alpha(colpal_DC[ind_ord_dis[inds_col][j]], .5), ...)
      graphics::points(
        (ysteps ~ xsteps),
        cex = .1,
        col = scales::alpha(colpal_DC[ind_ord_dis[inds_col][j]], .5)
      )

    }

  } else if (type=="boxplot"){

    graphics::boxplot(
      suby ~ subg,
      col = scales::alpha(colpal_DC[ind_ord_dis[inds_col]], .7),
      main = paste0("Boxplots colored by DC\n",main_title),
      ylab = "y",
      xlab = "Group index", ...
    )
  }else{
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par))

    graphics::par(mfrow=c(1,2))
    plot(suby ~ jitter(subg),
         col  = colpal_DC[subDC],
         xlab = "Group index",
         ylab = "y",
         main = paste0("Observations colored by DC\n",main_title), ...
    )
    plot(suby ~ jitter(subg),
         col  = colpal_OC[subOC],
         xlab = "Group index",
         ylab = "y",
         main = paste0("Observations colored by OC\n",main_title), ...
    )
  }
}
#' @name estimate_partition
#'
#' @param x The result of a call to \code{\link{estimate_partition}}.
#' @param DC_num An integer or a vector of integers indicating which distributional clusters to plot.
#' @param type What type of plot should be drawn. Available types are \code{"boxplot"}, \code{"ecdf"}, and \code{"scatter"}.
#' @param alt_palette Logical, the color palette to be used. Default is \code{R} base colors (\code{alt_palette = FALSE}).
#' @param ... Additional graphical parameters to be passed to the \code{plot} function.
#'
#' @importFrom graphics abline lines points boxplot
#' @importFrom stats median
#' @importFrom grDevices colorRampPalette
#' @importFrom scales alpha
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
plot.partition_vi <- function(x,
                       DC_num = NULL,
                       type = c("ecdf", "boxplot", "scatter"),
                       alt_palette = FALSE,
                       ...) {

  type   <- match.arg(type)
  if (is.null(DC_num)) {
    DC_num <- 1:max(x$dis_level)
  }
  if(max(DC_num) > max(x$dis_level)){
    stop(paste0("There are less estimated DCs than requested.\n",
                "Please provide a number for DC_num between 1 and ", max(x$dis_level) ))
  }


  max_OC <- max(x$obs_level$OC)
  max_DC <- max(x$dis_level)
  if(alt_palette){
    colpal_DC <- colorRampPalette(brewer.pal(8, "Dark2"))(max_DC)
    colpal_OC <- colorRampPalette(brewer.pal(8, "Dark2"))(max_OC)
  }else{
    colpal_DC <- 1:max_DC
    colpal_OC <- 1:max_OC
  }


  dix <- rank(tapply(x$obs_level$Y, x$obs_level$DC, stats::median))

  ind_ord_dis <- sapply(x$dis_level, function(g) dix[g])
  ind_ord_obs <- dix[x$obs_level$DC]


  inds_row <- which(ind_ord_obs %in% DC_num)
  inds_col <- which(ind_ord_dis %in% DC_num)

  suby <- (x$obs_level$Y[inds_row])
  subg <- (x$obs_level$G[inds_row])
  subDC <- x$obs_level$DC[inds_row]
  subOC <- x$obs_level$OC[inds_row]

  if(length(DC_num) == max(x$dis_level)){
    main_title <- paste0("All distributional clusters")
  }else if( length(DC_num) > 8 ){
    main_title <- NULL
  }else{
    main_title <- paste0("Distributional clusters # ", paste(DC_num, collapse = ", "))
  }



  if (type == "ecdf") {
    X <- unique(x$obs_level$G[inds_row])
    Nj <- table(x$obs_level$G)
    ysteps <- (0:(Nj[inds_col][1] - 1)) /
      (Nj[inds_col][1] - 1)
    xsteps <- sort(suby[subg == X[1]])


    plot(
      ysteps ~ xsteps,
      type = "l",
      xlim = c(min(suby), max(suby)),
      ylim = c(0, 1),
      col = scales::alpha(colpal_DC[ind_ord_dis[inds_col][1]], .5),
      xlab = "y",
      ylab = "eCDF",
      main = paste0("eCDFs colored by DC\n",main_title), ...
    )
    graphics::points(
      (ysteps ~ xsteps),
      cex = .1,
      col = scales::alpha(colpal_DC[ind_ord_dis[inds_col][1]], .5)
    )
    graphics::abline(h = c(0, 1),
                     col = "gray",
                     lty = 3)


    for (j in 2:length(X)) {
      ysteps <- (0:(Nj[X[j]] - 1)) / (Nj[X[j]] - 1)
      xsteps <- sort(suby[subg == X[j]])

      graphics::lines((ysteps ~ xsteps),
                      col = scales::alpha(colpal_DC[ind_ord_dis[inds_col][j]], .5), ...)
      graphics::points(
        (ysteps ~ xsteps),
        cex = .1,
        col = scales::alpha(colpal_DC[ind_ord_dis[inds_col][j]], .5)
      )

    }

  } else if (type=="boxplot"){

    graphics::boxplot(
      suby ~ subg,
      col = scales::alpha(colpal_DC[ind_ord_dis[inds_col]], .7),
      main = paste0("Boxplots colored by DC\n",main_title),
      ylab = "y",
      xlab = "group", ...
    )
  }else{
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par))

    graphics::par(mfrow=c(1,2))
    plot(suby ~ jitter(subg),
         col=colpal_DC[subDC],
         xlab = "Group index",
         ylab = "y",
         main = paste0("Observations colored by DC\n",main_title), ...
    )
    plot(suby ~ jitter(subg),
         col=colpal_OC[subOC],
         xlab = "Group index",
         ylab = "y",
         main = paste0("Observations colored by OC\n",main_title), ...
    )
  }
}

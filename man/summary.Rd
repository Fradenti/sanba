% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_SANBA_summary_estimate_partitions.R
\name{summary}
\alias{summary}
\alias{summary.SANvi}
\alias{summary.SANmcmc}
\alias{print.summary_mcmc}
\alias{print.summary_vi}
\alias{plot.summary_mcmc}
\alias{plot.summary_vi}
\title{Summarize the Estimated Observational and Distributional Partition}
\usage{
\method{summary}{SANvi}(object, ordered = TRUE, ...)

\method{summary}{SANmcmc}(object, ordered = TRUE, add_burnin = 0, ncores = 0, ...)

\method{print}{summary_mcmc}(x, ...)

\method{print}{summary_vi}(x, ...)

\method{plot}{summary_mcmc}(
  x,
  DC_num = NULL,
  type = c("ecdf", "boxplot", "scatter"),
  alt_palette = FALSE,
  ...
)

\method{plot}{summary_vi}(
  x,
  DC_num = NULL,
  type = c("ecdf", "boxplot", "scatter"),
  alt_palette = FALSE,
  ...
)
}
\arguments{
\item{object}{Object of class \code{SANmcmc} (usually, the result of a call to \code{\link{fit_fiSAN}},
\code{\link{fit_fSAN}}, or \code{\link{fit_CAM}} with \code{method = "MCMC"}) or \code{SANvi}
(the result of a call to \code{\link{fit_fiSAN}},\code{\link{fit_fSAN}}, or \code{\link{fit_CAM}} with \code{method = "VI"}).}

\item{ordered}{Logical, if \code{TRUE} (default), the function sorts the distributional cluster labels reflecting the
increasing values of medians of the data assigned to each DC.}

\item{...}{Additional graphical parameters to be passed to the \code{plot} function.}

\item{add_burnin}{Integer (default = 0). Number of observations to discard as additional burn-in (only for \code{SANmcmc} objects).}

\item{ncores}{A parameter to pass to the \code{salso::salso()} function (only for \code{SANmcmc} objects). The number of CPU cores to use for parallel computing; a value of zero indicates the use of all cores on the system.}

\item{x}{The result of a call to \code{\link{summary}}.}

\item{DC_num}{An integer or a vector of integers indicating which distributional clusters to plot.}

\item{type}{What type of plot should be drawn. Available types are \code{"boxplot"}, \code{"ecdf"}, and \code{"scatter"}.}

\item{alt_palette}{Logical, the color palette to be used. Default is \code{R} base colors (\code{alt_palette = FALSE}).}
}
\value{
A list of class \code{summary_vi} or \code{summary_mcmc} containing
\itemize{
  \item \code{obs_level}: a data frame containing the data values, their group indexes, and the observational and distributional clustering assignments for each observation.
  \item \code{dis_level}: a vector with the distributional clustering assignment for each unit.
}
}
\description{
Given the output of a \code{sanba} model-fitting function, estimate the observational and distributional partitions using \code{\link[salso:salso]{salso::salso()}} for MCMC, and the maximum a posteriori estimate for VI.
}
\examples{
set.seed(123)
y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
g <- c(rep(1:6, each = 10))
out <- fit_fSAN(y = y, group = g, "VI", vi_param = list(n_runs = 10))
plot(out)
clust <- summary(out)
clust
plot(clust, lwd = 2, alt_palette = TRUE)
plot(clust, type = "scatter", alt_palette = FALSE, cex = 2)


set.seed(123)
y <- c(rnorm(40,0,0.3), rnorm(20,5,0.3))
g <- c(rep(1:6, each = 10))
out <- fit_fSAN(y = y, group = g, "MCMC", mcmc_param=list(nrep=500,burn=200))
plot(out)
clust <- summary(out)
clust
plot(clust, lwd = 2)
plot(clust,  type = "boxplot", alt_palette = TRUE)
plot(clust,  type = "scatter", alt_palette = TRUE, cex = 2, pch = 4)

}
\seealso{
\code{\link[salso:salso]{salso::salso()}}, \code{\link{print.SANmcmc}}, \code{\link{plot.SANmcmc}}
}

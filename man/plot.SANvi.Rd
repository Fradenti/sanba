% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_SANBA_plotting.R
\name{plot.SANvi}
\alias{plot.SANvi}
\title{Visual Check of the Convergence of the VI Output}
\usage{
\method{plot}{SANvi}(x, ...)
}
\arguments{
\item{x}{Object of class \code{SANvi} (usually, the result of a call to \code{fit_CAM}, \code{fit_fiSAN}, or \code{fit_fSAN}, used with the \code{est_method = "VI"} argument).}

\item{...}{Ignored.}
}
\value{
The function plots the path followed by the ELBO and its subsequent differences.
}
\description{
Plot method for objects of class \code{SANvi}.
The function displays two graphs. The left plot shows the progression of all the ELBO values as a function of the iterations.
The right plots shows the ELBO increments between successive iterations of the best run on a log scale (note: increments should always be positive).
}
\examples{
set.seed(123)
y <- c(rnorm(200,0,0.3), rnorm(100,5,0.3))
g <- c(rep(1,150), rep(2, 150))
out <- fit_fSAN(y = y, group = g, "VI", vi_param = list(n_runs = 2))
plot(out)
}

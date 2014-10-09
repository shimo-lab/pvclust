\name{msfit}
\alias{msfit}
\alias{plot.msfit}
\alias{lines.msfit}
\alias{summary.msfit}
\title{Curve Fitting for Multiscale Bootstrap Resampling}
\description{\code{msfit} performs curve fitting for multiscale
  bootstrap resampling. It generates an object of class
  \code{msfit}. Several generic methods are available.
}
\usage{
msfit(bp, r, nboot)

\method{plot}{msfit}(x, curve=TRUE, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, ...)

\method{lines}{msfit}(x, col=2, lty=1, ...)

\method{summary}{msfit}(object, digits=3, ...)
}
\arguments{
  \item{bp}{numeric vector of bootstrap probability values.}
  \item{r}{numeric vector of relative sample size of bootstrap samples
    defined as \eqn{r=n'/n} for original sample size \eqn{n} and
    bootstrap sample size \eqn{n'}.}
  \item{nboot}{numeric value (vector) of the number of bootstrap
    replications.}
  \item{x}{object of class \code{msfit}.}
  \item{curve}{logical. If \code{TRUE}, the fitted curve is drawn.}
  \item{main, sub, xlab, ylab, col, lty}{generic graphic parameters.}
  \item{object}{object of class \code{msfit}.}
  \item{digits}{integer indicating the precision to be used in rounding.}
  \item{...}{other parameters to be used in the functions.}
}
\details{
  function \code{msfit} performs the curve fitting for multiscale
  bootstrap resampling. In package \code{pvclust} this function is only
  called from the function \code{pvclust} (or \code{parPvclust}), and
  may never be called from users. However one can access a list of
  \code{msfit} objects by \code{x$msfit}, where \code{x} is an object of
  class \code{pvclust}.
}
\value{\code{msfit} returns an object of class \code{msfit}. It contains
  the following objects:
  \item{p}{numeric vector of \eqn{p}-values. \code{au} is AU
    (Approximately Unbiased) \eqn{p}-value computed by multiscale
    bootstrap resampling, which is more accurate than BP value
    (explained below) as unbiased \eqn{p}-value. \code{bp} is BP
    (Bootstrap Probability) value, which is simple but tends to be
    unbiased when the absolute value of \code{c} (a value in \code{coef}
    vector, explained below) is large.}
  \item{se}{numeric vector of estimated standard errors of \eqn{p}-values.}
  \item{coef}{numeric vector related to geometric aspects of
    hypotheses. \code{v} is signed distance and \code{c} is curvature of
    the boundary.}
  \item{df}{numeric value of the degree of freedom in curve fitting.}
  \item{rss}{residual sum of squares.}
  \item{pchi}{\eqn{p}-value of chi-square test based on asymptotic theory.}
}
\references{
  Shimodaira, H. (2004)
  "Approximately unbiased tests of regions using multistep-multiscale
  bootstrap resampling",
  \emph{Annals of Statistics}, 32, 2616-2641.

  Shimodaira, H. (2002)
  "An approximately unbiased test of phylogenetic tree selection",
  \emph{Systematic Biology}, 51, 492-508.
}
\author{Ryota Suzuki \email{suzuki@ef-prime.com}}
\keyword{htest}
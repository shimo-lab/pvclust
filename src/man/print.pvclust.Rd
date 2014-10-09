\name{print.pvclust}
\alias{print.pvclust}
\title{Print Function for Pvclust Object}
\description{print clustering method and distance measure used in
  hierarchical clustering, \eqn{p}-values and related statistics for
  a \code{pvclust} object.
}
\usage{
\method{print}{pvclust}(x, which=NULL, digits=3, ...)
}
\arguments{
  \item{x}{object of class \code{pvclust}.}
  \item{which}{numeric vector which specifies the numbers of edges
    (clusters) of which the values are printed. If \code{NULL} is given,
    it prints the values of all edges. The default is \code{NULL}.}
  \item{digits}{integer indicating the precision to be used in
    rounding.}
  \item{...}{other parameters used in the function.}
}
\value{this function prints \eqn{p}-values and some related
  statistics.
  \item{au}{AU (Approximately Unbiased) \eqn{p}-value, which is more
    accurate than BP value as unbiased \eqn{p}-value. It is computed by
    multiscale bootstrap resampling.}
  \item{bp}{BP (Bootstrap Probability) value, which is a simple
    statistic computed by bootstrap resampling. This value tends to be
    biased as \eqn{p}-value when the absolute value of \code{c}
    (explained below) is large.}
  \item{se.au, se.bp}{estimated standard errors for \code{au} and
    \code{bp}, respectively.}
  \item{v, c}{values related to geometric aspects of
    hypotheses. \code{v} is signed distance and \code{c} is curvature of
    the boundary.}
  \item{pchi}{\eqn{p}-values of chi-square test based on asymptotic theory.}
}
\author{Ryota Suzuki \email{suzuki@ef-prime.com}}
\keyword{print}
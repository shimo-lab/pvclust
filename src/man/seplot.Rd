\name{seplot}
\alias{seplot}
\title{Diagnostic Plot for Standard Error of p-value}
\description{draws diagnostic plot for standard error of \eqn{p}-value
  for pvclust object.
}
\usage{
seplot(object, type=c("au", "bp"), identify=FALSE, main=NULL,
xlab=NULL, ylab=NULL, ...)
}
\arguments{
  \item{object}{object of class \code{pvclust}.}
  \item{type}{the type of \eqn{p}-value to be plotted, one of \code{"au"}
    or \code{"bp"}.}
  \item{identify}{logical. If \code{TRUE}, edge numbers can be identified
    interactively. See \code{\link{identify}} for basic usage.}
  \item{main, xlab, ylab}{generic graphic parameters. See
    \code{\link[graphics]{par}} for details.}
  \item{...}{other graphical parameters to be passed to generic
    \code{plot} or \code{identify} function.}
}
\author{Ryota Suzuki \email{suzuki@ef-prime.com}}
\keyword{hplot}

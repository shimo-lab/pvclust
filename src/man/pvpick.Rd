\name{pvpick}
\alias{pvpick}
\alias{pvrect}
\alias{lines.pvclust}
\title{Find Clusters with High/Low P-values}
\description{find clusters with relatively high/low
  \eqn{p}-values. \code{pvrect} and \code{lines} (S3 method for class
  \code{pvclust}) highlight such clusters in existing plot, and
  \code{pvpick} returns a list of such clusters.
}
\usage{
pvpick(x, alpha=0.95, pv="au", type="geq", max.only=TRUE)

pvrect(x, alpha=0.95, pv="au", type="geq", max.only=TRUE, border=2, ...)

\method{lines}{pvclust}(x, alpha=0.95, pv="au", type="geq", col=2, lwd=2, ...)
}
\arguments{
  \item{x}{object of class \code{pvclust}.}
  \item{alpha}{threshold value for \eqn{p}-values.}
  \item{pv}{character string which specifies the \eqn{p}-value to be
    used. It should be either of \code{"au"} or \code{"bp"},
    corresponding to AU \eqn{p}-value or BP value, respectively. See
    \code{plot.pvclust} for details.}
  \item{type}{one of \code{"geq"}, \code{"leq"}, \code{"gt"} or
    \code{"lt"}. If \code{"geq"} is specified, clusters with
    \eqn{p}-value \emph{greater than or equals} the threshold given by
    \code{"alpha"} are returned or displayed. Likewise \code{"leq"}
    stands for \emph{lower than or equals}, \code{"gt"} for
    \emph{greater than} and \code{"lt"} for \emph{lower than} the
    threshold value. The default is \code{"geq"}.}
  \item{max.only}{logical. If some of clusters with high/low
    \eqn{p}-values have inclusion relation, only the largest cluster is
    returned (or displayed) when \code{max.only=TRUE}.}
  \item{border}{numeric value which specifies the color of borders of
    rectangles.}
  \item{col}{numeric value which specifies the color of lines.}
  \item{lwd}{numeric value which specifies the width of lines.}
  \item{...}{other graphic parameters to be used.}
}
\value{\code{pvpick} returns a list which contains the following values.
  \item{clusters}{a list of character string vectors. Each vector
    corresponds to the names of objects in each cluster.}
  \item{edges}{numeric vector of edge numbers. The \eqn{i}'th element
    (number) corresponds to the \eqn{i}'th name vector in
    \code{clusters}.}
  }
\author{Ryota Suzuki \email{suzuki@ef-prime.com}}
\keyword{aplot}
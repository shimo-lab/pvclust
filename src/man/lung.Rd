\name{lung}
\docType{data}
\alias{lung}
\title{DNA Microarray Data of Lung Tumors}
\description{
DNA Microarray data of 73 lung tissues including 67 lung tumors.
There are 916 observations of genes for each lung tissue.
}
\usage{data(lung)}
\format{data frame of size \eqn{916 \times 73}{(916, 73)}.}
\details{
  This dataset has been modified from original data. Each one
  observation of duplicate genes has been removed. See \code{source}
  section in this help for original data source.
}
\examples{
\donttest{
## Reading the data
data(lung)

## Multiscale Bootstrap Resampling
lung.pv <- pvclust(lung, nboot=100)

## CAUTION: nboot=100 may be too small for actual use.
##          We suggest nboot=1000 or larger.
##          plot/print functions will be useful for diagnostics.

## Plot the result
plot(lung.pv, cex=0.8, cex.pv=0.7)

ask.bak <- par()$ask
par(ask=TRUE)

pvrect(lung.pv, alpha=0.9)
msplot(lung.pv, edges=c(51,62,68,71))

par(ask=ask.bak)

## Print a cluster with high p-value
lung.pp <- pvpick(lung.pv, alpha=0.9)
lung.pp$clusters[[2]]

## Print its edge number
lung.pp$edges[2]
}

## We recommend parallel computing for large dataset as this one
\dontrun{
library(snow)
cl <- makeCluster(10, type="MPI")
lung.pv <- parPvclust(cl, lung, nboot=1000)
}
}
\source{\url{http://genome-www.stanford.edu/lung_cancer/adeno/}}
\references{
  Garber, M. E. et al. (2001)
  "Diversity of gene expression in adenocarcinoma of the lung",
  \emph{Proceedings of the National Academy of Sciences},
  98, 13784-13789.
}
\keyword{datasets}

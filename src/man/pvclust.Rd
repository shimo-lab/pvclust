\name{pvclust}
\alias{pvclust}
\alias{parPvclust}
\title{Calculating P-values for Hierchical Clustering}
\description{
  calculates \eqn{p}-values for hierarchical clustering via
  multiscale bootstrap resampling. Hierarchical clustering is done for
  given data and \eqn{p}-values are computed for each of the clusters.
}
\usage{
pvclust(data, method.hclust="average",
        method.dist="correlation", use.cor="pairwise.complete.obs",
        nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE)

parPvclust(cl=NULL, data, method.hclust="average",
           method.dist="correlation", use.cor="pairwise.complete.obs",
           nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE,
           init.rand=TRUE, seed=NULL, iseed=NULL)

}
\arguments{
  \item{data}{numeric data matrix or data frame.}
  \item{method.hclust}{
    the agglomerative method used in hierarchical clustering. This
    should be (an abbreviation of) one of \code{"average"}, \code{"ward.D"},
    \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"mcquitty"},
    \code{"median"} or \code{"centroid"}. The default is
    \code{"average"}. See \code{method} argument in
    \code{\link[stats]{hclust}}.
  }
  \item{method.dist}{the distance measure to be used. This should be
    a character string, or a function which returns a \code{dist} object.
    A character string should be (an abbreviation of)
    one of \code{"correlation"}, \code{"uncentered"},
    \code{"abscor"} or those which are allowed for \code{method}
    argument in \code{dist} function. The default is
    \code{"correlation"}. See \emph{details} section in this help and
    \code{method} argument in \code{\link[stats]{dist}}.
    }
  \item{use.cor}{character string which specifies the method for
    computing correlation with data including missing values. This
    should be (an abbreviation of) one of \code{"all.obs"},
    \code{"complete.obs"} or \code{"pairwise.complete.obs"}. See
    the \code{use} argument in \code{\link[stats]{cor}} function.
  }
  \item{nboot}{the number of bootstrap replications. The default is
    \code{1000}.}
  \item{r}{numeric vector which specifies the relative sample sizes of
    bootstrap replications. For original sample size \eqn{n} and
    bootstrap sample size \eqn{n'}, this is defined as \eqn{r=n'/n}.}
  \item{store}{locical. If \code{store=TRUE}, all bootstrap replications
  are stored in the output object. The default is \code{FALSE}.}
  \item{cl}{a cluster object created by package \pkg{parallel} or \pkg{snow}.
    If NULL, use the registered default cluster.}
  \item{weight}{logical. If \code{weight=TRUE}, resampling is made by
    weight vector instead of index vector. Useful for large \code{r}
    value (\code{r>10}).  Currently, available only for distance
    \code{"correlation"} and \code{"abscor"}.}
%  \item{init.rand}{logical. If \code{init.rand=TRUE}, random number
%    generators are initialized at child processes. Random seeds can be
%    set by \code{seed} argument.}
%  \item{seed}{integer vector of random seeds. It should have the same
%    length as \code{cl}. If \code{NULL} is specified,
%    \code{1:length(cl)} is used as seed vector. The default is \code{NULL}.}
  \item{init.rand}{logical. If \code{init.rand=TRUE}, random number generators are initialized.
    Use \code{iseed} argument to achieve reproducible results.}
  \item{seed}{integer vector of random seeds. It should have the same
    length as \code{cl}. \strong{This argument is duplicated and will be unavailable in the future. Consider using \code{iseed} instead.} }
  \item{iseed}{an integer to be passed to \code{clusterSetRNGStream} when \code{init.rand=TRUE}.}

}
\details{
  Function \code{pvclust} conducts multiscale bootstrap resampling to calculate
  \eqn{p}-values for each cluster in the result of hierarchical
  clustering. \code{parPvclust} is the parallel version of this
  procedure which depends on package \pkg{parallel} for parallel computation.

  For data expressed as \eqn{(n \times p)}{(n, p)} matrix or data frame, we
  assume that the data is \eqn{n} observations of \eqn{p} objects, which
  are to be clustered. The \eqn{i}'th row vector corresponds to the
  \eqn{i}'th observation of these objects and the \eqn{j}'th column
  vector corresponds to a sample of \eqn{j}'th object with size \eqn{n}.

  There are several methods to measure the dissimilarities between
  objects. For data matrix \eqn{X=\{x_{ij}\}}{X},
  \code{"correlation"}
  method takes
  \deqn{
    1 - \frac{
      \sum_{i=1}^n (x_{ij} - \bar{x}_j) (x_{ik} - \bar{x}_k)
    }
    {
      \sqrt{\sum_{i=1}^n (x_{ij} - \bar{x}_j)^2}
      \sqrt{\sum_{i=1}^n (x_{ik} - \bar{x}_k)^2}
    }
    }{%
      1 - cor(X)[j,k]
    }
  for dissimilarity between \eqn{j}'th and \eqn{k}'th object, where
  \eqn{\bar{x}_j = \frac{1}{n} \sum_{i=1}^n x_{ij} \mbox{and}
  \bar{x}_k = \frac{1}{n} \sum_{i=1}^n x_{ik}}{cor is function \code{cor}}.

  \code{"uncentered"} takes uncentered sample correlation
  \deqn{
    1 - \frac{
      \sum_{i=1}^n x_{ij} x_{ik}
    }
    {
      \sqrt{\sum_{i=1}^n x_{ij}^2}
      \sqrt{\sum_{i=1}^n x_{ik}^2}
    }
    }{%
      1 - sum(x[,j] * x[,k]) / (sqrt(sum(x[,j]^2)) * sqrt(sum(x[,k]^2)))
    }
  and \code{"abscor"} takes the absolute value of sample correlation
  \deqn{
    1 - \ \Biggl| \frac{
      \sum_{i=1}^n (x_{ij} - \bar{x}_j) (x_{ik} - \bar{x}_k)
    }
    {
      \sqrt{\sum_{i=1}^n (x_{ij} - \bar{x}_j)^2}
      \sqrt{\sum_{i=1}^n (x_{ik} - \bar{x}_k)^2}
    } \Biggl|.
  }{%
    1 - abs(cor(X)[j,k]).
  }
  }
\value{
  \item{hclust}{hierarchical clustering for original data generated by
    function \code{hclust}. See \code{\link[stats]{hclust}} for details.}
  \item{edges}{data frame object which contains \eqn{p}-values and
    supporting informations such as standard errors.}
  \item{count}{data frame object which contains primitive information
    about the result of multiscale bootstrap resampling.}
  \item{msfit}{list whose elements are results of curve fitting for
    multiscale bootstrap resampling, of class \code{msfit}. See
    \code{\link{msfit}} for details.}
  \item{nboot}{numeric vector of number of bootstrap replications.}
  \item{r}{numeric vector of the relative sample size for bootstrap
    replications.}
  \item{store}{list contains bootstrap replications if \code{store=TRUE}
  was given for function \code{pvclust} or \code{parPvclust}.}
}
\seealso{\code{\link{lines.pvclust}}, \code{\link{print.pvclust}},
  \code{\link{msfit}}, \code{\link{plot.pvclust}},
  \code{\link{text.pvclust}}, \code{\link{pvrect}} and
  \code{\link{pvpick}}.}
\references{
  Shimodaira, H. (2004)
  "Approximately unbiased tests of regions using multistep-multiscale
  bootstrap resampling",
  \emph{Annals of Statistics}, 32, 2616-2641.

  Shimodaira, H. (2002)
  "An approximately unbiased test of phylogenetic tree selection",
  \emph{Systematic Biology}, 51, 492-508.

  Suzuki, R. and Shimodaira, H. (2004)
  "An application of multiscale bootstrap resampling to hierarchical
  clustering of microarray data: How accurate are these clusters?",
  \emph{The Fifteenth International Conference on Genome Informatics 2004},
  P034.

  \url{http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/pvclust/}
}
\examples{
### example using Boston data in package MASS
data(Boston, package = "MASS")

## multiscale bootstrap resampling
boston.pv <- pvclust(Boston, nboot=100)

## CAUTION: nboot=100 may be too small for actual use.
##          We suggest nboot=1000 or larger.
##          plot/print functions will be useful for diagnostics.

## plot dendrogram with p-values
plot(boston.pv)

ask.bak <- par()$ask
par(ask=TRUE)

## highlight clusters with high au p-values
pvrect(boston.pv)

## print the result of multiscale bootstrap resampling
print(boston.pv, digits=3)

## plot diagnostic for curve fitting
msplot(boston.pv, edges=c(2,4,6,7))

par(ask=ask.bak)

## print clusters with high p-values
boston.pp <- pvpick(boston.pv)
boston.pp

### Using a custom distance measure

## Define a distance function which returns an object of class "dist".
## The function must have only one argument "x" (data matrix or data.frame).
cosine <- function(x) {
    x <- as.matrix(x)
    y <- t(x) \%*\% x
    res <- 1 - y / (sqrt(diag(y)) \%*\% t(sqrt(diag(y))))
    res <- as.dist(res)
    attr(res, "method") <- "cosine"
    return(res)
}

result <- pvclust(Boston, method.dist=cosine, nboot=100)
plot(result)

\dontrun{
### parallel computation
library(parallel)
cl <- makeCluster(2, type = "PSOCK")

## parallel version of pvclust
boston.pv <- parPvclust(cl, Boston, nboot=1000)
stopCluster(cl)
}
}
\author{Ryota Suzuki \email{suzuki@ef-prime.com}}
\keyword{cluster}

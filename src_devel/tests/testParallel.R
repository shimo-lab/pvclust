library(pvclust)
library(MASS)

pv <- pvclust(Boston, nboot = 1e+4)

## parallel version
library(parallel)
n <- max(2, detectCores())
cl <- makeCluster(n, type="PSOCK")
pv.par <- parPvclust(cl, Boston, nboot=1e+4)
stopCluster(cl)

compare.pvclust <- function(x, y, bp.threshold = 0.01, au.threshold = 0.01) {

  ## check if both have the same hclust object
  stopifnot(identical(x$hclust, y$hclust))
  
  e.x <- x$edges
  e.y <- y$edges
  
  stopifnot(nrow(e.x) == nrow(e.y))
  
  bp.error <- abs(e.x$bp - e.y$bp)
  au.error <- abs(e.x$au - e.y$au)
  
  stopifnot(mean(bp.error) <= bp.threshold)
  if(max(bp.error) > bp.threshold) warning("max(bp.error) = ", max(bp.error))
  
  stopifnot(mean(au.error) <= au.threshold)
  if(max(au.error) > au.threshold) warning("max(au.error) = ", max(au.error))
  
}

compare.pvclust(x = pv, y = pv.par)
### internal function for non-parallel pvclust
pvclust.nonparallel <- function(data, method.hclust, method.dist, use.cor, nboot, r,
                                store, weight, iseed, quiet)
{
  # initialize random seed
  if(!is.null(iseed))
    set.seed(seed = iseed)
  
  # data: (n,p) matrix, n-samples, p-variables
  n <- nrow(data); p <- ncol(data)
  
  # hclust for original data
  #    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
  #                 "median", "centroid")
  #    method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
  
  # Use custom distance function
  if(is.function(method.dist)) {
    distance <- method.dist(data)
  } else {
    distance <- dist.pvclust(data, method=method.dist, use.cor=use.cor)
  }
  
  data.hclust <- hclust(distance, method=method.hclust)
  
  # ward -> ward.D
  # only if R >= 3.1.0
  if(method.hclust == "ward" && getRversion() >= '3.1.0') {
      method.hclust <- "ward.D"
  }

  # multiscale bootstrap
  size <- floor(n*r)
  rl <- length(size)
  
  if(rl == 1) {
    if(r != 1.0)
      warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
    
    r <- list(1.0)
  }
  else
    r <- as.list(size/n)
  
  mboot <- lapply(r, boot.hclust, data=data, object.hclust=data.hclust, nboot=nboot,
                  method.dist=method.dist, use.cor=use.cor,
                  method.hclust=method.hclust, store=store, weight=weight, quiet=quiet)
  
  result <- pvclust.merge(data=data, object.hclust=data.hclust, mboot=mboot)
  
  return(result)
}


### internal function for parallel pvclust
pvclust.parallel <- function(cl, data, method.hclust, method.dist, use.cor,
                             nboot, r, store, weight, init.rand=NULL, iseed, quiet,
                             parallel.check)
{
  if(parallel.check) {    
    check.result <- check.parallel(cl=cl, nboot=nboot)
    if(!check.result) {
      msg <- paste(attr(check.result, "msg"), ". non-parallel version is executed", sep = "")
      warning(msg)
      return(pvclust.nonparallel(data=data, method.hclust=method.hclust, method.dist=method.dist,
                                 use.cor=use.cor, nboot=nboot, r=r, store=store, weight=weight, iseed=iseed, quiet=quiet))
    }
  }
  
  # check package versions
  pkg.ver <-parallel::clusterCall(cl, packageVersion, pkg = "pvclust")
  r.ver   <- parallel::clusterCall(cl, getRversion)
  if(length(unique(pkg.ver)) > 1 || length(unique(r.ver)) > 1) {
    
    node.name <- parallel::clusterEvalQ(cl, Sys.info()["nodename"])
    version.table <- data.frame(
      node=seq_len(length(node.name)),
      name=unlist(node.name),
      R=unlist(lapply(r.ver, as.character)),
      pvclust=unlist(lapply(pkg.ver, as.character)))
    
    if(nrow(version.table) > 10)
      table.out <- c(capture.output(print(head(version.table, n=10), row.names=FALSE)), "    ...")
    else
      table.out <- capture.output(print(version.table, row.names=FALSE))
    
    warning("R/pvclust versions are not unique:\n",
      paste(table.out, collapse="\n"))
  }
  
  if(!is.null(init.rand))
    warning("\"init.rand\" option is deprecated. It is available for back compatibility but will be unavailable in the future.\nSpecify a non-NULL value of \"iseed\" to initialize random seed.")
  
  #   if(init.rand) {
  #     if(is.null(iseed) && !is.null(seed)) {
  #       warning("\"seed\" option is deprecated. It is available for back compatibility but will be unavailable in the future.\nConsider using \"iseed\" instead.")
  #       
  #       if(length(seed) != length(cl))
  #         stop("seed and cl should have the same length.")
  #       
  #       # setting random seeds
  #       parallel::parLapply(cl, as.list(seed), set.seed)
  #     } else {
  #       parallel::clusterSetRNGStream(cl = cl, iseed = iseed)
  #     }
  #   }
  
  if(!is.null(iseed) && (is.null(init.rand) || init.rand))
    parallel::clusterSetRNGStream(cl = cl, iseed = iseed)
  
  # data: (n,p) matrix, n-samples, p-variables
  n <- nrow(data); p <- ncol(data)
  
  # hclust for original data
  if(is.function(method.dist)) {
    # Use custom distance function
    distance <- method.dist(data)
  } else {
    distance <- dist.pvclust(data, method=method.dist, use.cor=use.cor)
  }
  
  data.hclust <- hclust(distance, method=method.hclust)
  
  # ward -> ward.D
  # only if R >= 3.1.0
  if(method.hclust == "ward" && getRversion() >= '3.1.0') {
    method.hclust <- "ward.D"
  }
  
  # multiscale bootstrap
  size <- floor(n*r)
  rl <- length(size)
  
  if(rl == 1) {
    if(r != 1.0)
      warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
    
    r <- list(1.0)
  }
  else
    r <- as.list(size/n)
  
  ncl <- length(cl)
  nbl <- as.list(rep(nboot %/% ncl,times=ncl))
  
  if((rem <- nboot %% ncl) > 0)
    nbl[1:rem] <- lapply(nbl[1:rem], "+", 1)
  
  if(!quiet)
    cat("Multiscale bootstrap... ")
  
  mlist <- parallel::parLapply(cl, nbl, pvclust.node,
                               r=r, data=data, object.hclust=data.hclust, method.dist=method.dist,
                               use.cor=use.cor, method.hclust=method.hclust,
                               store=store, weight=weight, quiet=quiet)
  if(!quiet)
    cat("Done.\n")
  
  mboot <- mlist[[1]]
  
  for(i in 2:ncl) {
    for(j in 1:rl) {
      mboot[[j]]$edges.cnt <- mboot[[j]]$edges.cnt + mlist[[i]][[j]]$edges.cnt
      mboot[[j]]$nboot <- mboot[[j]]$nboot + mlist[[i]][[j]]$nboot
      mboot[[j]]$store <- c(mboot[[j]]$store, mlist[[i]][[j]]$store)
    }
  }
  
  result <- pvclust.merge( data=data, object.hclust=data.hclust, mboot=mboot)
  
  return(result)
}

hc2axes <- function(x)
{
  A <- x$merge # (n,n-1) matrix
  n <- nrow(A) + 1
  x.axis <- c()
  y.axis <- x$height
  
  x.tmp  <- rep(0,2)
  zz     <- match(1:length(x$order),x$order)
  
  for(i in 1:(n-1)) {
    ai <- A[i,1]
    
    if(ai < 0)
      x.tmp[1] <- zz[-ai]
    else
      x.tmp[1] <- x.axis[ai]
    
    ai <- A[i,2]
    
    if(ai < 0)
      x.tmp[2] <- zz[-ai]
    else
      x.tmp[2] <- x.axis[ai]
    
    x.axis[i] <- mean(x.tmp)
  }
  
  return(data.frame(x.axis=x.axis,y.axis=y.axis))
}

hc2split <- function(x)
{
  A <- x$merge # (n-1,n) matrix
  n <- nrow(A) + 1
  B <- list()
  
  for(i in 1:(n-1)){
    ai <- A[i,1]
    
    if(ai < 0)
      B[[i]] <- -ai
    else
      B[[i]] <- B[[ai]]        
    
    ai <- A[i,2]
    
    if(ai < 0)
      B[[i]] <- sort(c(B[[i]],-ai))
    else
      B[[i]] <- sort(c(B[[i]],B[[ai]]))
  }
  
  CC <- matrix(rep(0,n*(n-1)),nrow=(n-1),ncol=n)
  
  for(i in 1:(n-1)){
    bi <- B[[i]]
    m <- length(bi)
    for(j in 1:m)
      CC[i,bi[j]] <- 1
  }
  
  split <- list(pattern=apply(CC,1,paste,collapse=""), member=B)
  
  return(split)
}

pvclust.node <- function(x, r, ...)
{
  #    require(pvclust)
  mboot.node <- lapply(r, boot.hclust, nboot=x, ...)
  return(mboot.node)
}

boot.hclust <- function(r, data, object.hclust, method.dist, use.cor,
                        method.hclust, nboot, store, weight=FALSE, quiet=FALSE)
{ 
  n     <- nrow(data)
  size  <- round(n*r, digits=0)
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n
  
  pattern   <- hc2split(object.hclust)$pattern
  edges.cnt <- table(factor(pattern)) - table(factor(pattern))
  st <- list()
  
  # bootstrap start
  rp <- as.character(round(r,digits=2)); if(r == 1) rp <- paste(rp,".0",sep="")
  if(!quiet)
      cat(paste("Bootstrap (r = ", rp, ")... ", sep=""))
  w0 <- rep(1,n) # equal weight
  na.flag <- 0
  
  for(i in 1:nboot){
    if(weight && r>10) {  ## <- this part should be improved
      w1 <- as.vector(rmultinom(1,size,w0)) # resampled weight
      suppressWarnings(distance <- distw.pvclust(data,w1,method=method.dist,use.cor=use.cor))
    } else {
      smpl <- sample(1:n, size, replace=TRUE)
      if(is.function(method.dist)) {
        suppressWarnings(distance  <- method.dist(data[smpl,]))
      } else {
        suppressWarnings(distance  <- dist.pvclust(data[smpl,],method=method.dist,use.cor=use.cor))
      }
    }
    if(all(is.finite(distance))) { # check if distance is valid
      x.hclust  <- hclust(distance,method=method.hclust)
      pattern.i <- hc2split(x.hclust)$pattern # split
      edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
    } else {
      x.hclust <- NULL
      na.flag <- 1
    }
    
    if(store)
      st[[i]] <- x.hclust
  }
  if(!quiet)
    cat("Done.\n")
  # bootstrap done
  
  if(na.flag == 1)
    warning(paste("inappropriate distance matrices are omitted in computation: r = ", r), call.=FALSE)
  
  boot <- list(edges.cnt=edges.cnt, method.dist=method.dist, use.cor=use.cor,
               method.hclust=method.hclust, nboot=nboot, size=size, r=r, store=st)
  class(boot) <- "boot.hclust"
  
  return(boot)
}

pvclust.merge <- function(data, object.hclust, mboot){
  
  pattern <- hc2split(object.hclust)$pattern
  
  r     <- unlist(lapply(mboot,"[[","r"))
  nboot <- unlist(lapply(mboot,"[[","nboot"))
  store <- lapply(mboot,"[[", "store")
  
  rl <- length(mboot)
  ne <- length(pattern)
  
  edges.bp <- edges.cnt <- data.frame(matrix(rep(0,ne*rl),nrow=ne,ncol=rl))
  row.names(edges.bp) <- pattern
  names(edges.cnt) <- paste("r", 1:rl, sep="")
  
  for(j in 1:rl) {
    edges.cnt[,j] <- as.vector(mboot[[j]]$edges.cnt) 
    edges.bp[,j]  <- edges.cnt[,j] / nboot[j]
  }
  
  ms.fitted <- lapply(as.list(1:ne),
                      function(x, edges.bp, r, nboot){
                        msfit(as.vector(t(edges.bp[x,])), r, nboot)},
                      edges.bp, r, nboot)
  class(ms.fitted) <- "mslist"
  
  p    <- lapply(ms.fitted,"[[","p")
  se   <- lapply(ms.fitted,"[[","se")
  coef <- lapply(ms.fitted,"[[","coef")
  
  au    <- unlist(lapply(p,"[[","au"))
  bp    <- unlist(lapply(p,"[[","bp"))
  se.au <- unlist(lapply(se,"[[","au"))
  se.bp <- unlist(lapply(se,"[[","bp"))
  v     <- unlist(lapply(coef,"[[","v"))
  cc    <- unlist(lapply(coef,"[[","c"))
  pchi  <- unlist(lapply(ms.fitted,"[[","pchi"))
  
  edges.pv <- data.frame(au=au, bp=bp, se.au=se.au, se.bp=se.bp,
                         v=v, c=cc, pchi=pchi)
  
  row.names(edges.pv) <- row.names(edges.cnt) <- 1:ne
  
  result <- list(hclust=object.hclust, edges=edges.pv, count=edges.cnt,
                 msfit=ms.fitted, nboot=nboot, r=r, store=store)
  
  class(result) <- "pvclust"
  return(result)
}

dist.pvclust <- function(x, method="euclidean", use.cor="pairwise.complete.obs")
{
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - cor(x, method="pearson", use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(cor(x,method="pearson",use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  else if(!is.na(pmatch(method,"uncentered"))){
    if(sum(is.na(x)) > 0){
      x <- na.omit(x)
      warning("Rows including NAs were omitted")
    }
    x  <- as.matrix(x)
    P  <- crossprod(x)
    qq <- matrix(diag(P),ncol=ncol(P))
    Q  <- sqrt(crossprod(qq))
    res <- as.dist(1 - P/Q)
    attr(res,"method") <- "uncentered"
    return(res)
  }
  else
    dist(t(x),method)
}


corw <- function(x,w,
                 use=c("all.obs","complete.obs","pairwise.complete.obs")
) {
  if(is.data.frame(x)) x <- as.matrix(x)
  x <- x[w>0,,drop=F]
  w <- w[w>0]
  
  n <- nrow(x) # sample size
  m <- ncol(x) # number of variables
  if(missing(w)) w <- rep(1,n)
  r <- matrix(0,m,m,dimnames=list(colnames(x),colnames(x)))
  diag(r) <- 1
  use <- match.arg(use)
  
  pairu <- F
  if(use=="all.obs") {
    u <- rep(T,n)
  } else if(use=="complete.obs") {
    u <- apply(x,1,function(y) !any(is.na(y)))
  } else if(use=="pairwise.complete.obs") {
    pairu <- T
    ux <- is.finite(x)
  } else stop("unknown use")
  
  for(i in 1+seq(length=m-1)) {
    for(j in seq(length=i-1)) {
      if(pairu) u <- ux[,i] & ux[,j]
      wu <- w[u]; xi <- x[u,i]; xj <- x[u,j]
      ws <- sum(wu)
      if(ws > 1e-8) {
        xi <- xi - sum(wu*xi)/ws
        xj <- xj - sum(wu*xj)/ws
        vxi <- sum(wu*xi*xi)/ws
        vxj <- sum(wu*xj*xj)/ws
        if(min(vxi,vxj) > 1e-8)  {
          vxij <- sum(wu*xi*xj)/ws
          rij <- vxij/sqrt(vxi*vxj)
        } else {
          rij <- 0
        }
      } else {
        rij <- 0
      }
      r[i,j] <- r[j,i] <- rij
    }
  }
  r
}

### calculate distance by weight
distw.pvclust <- function(x,w,method="correlation", use.cor="pairwise.complete.obs")
{
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - corw(x,w, use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(corw(x,w, use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  stop("wrong method")
}

### check whether parallel computation is appropriate
check.parallel <- function(cl, nboot) {
  res <- FALSE
  
### will be used when defaultCluster(cl) becomes publicly available
#   # check whether cl is a cluster, or a default cluster is available
#   if(!inherits(cl, "cluster")) {
#     try_result <- try(cl <- parallel:::defaultCluster(cl), silent=TRUE)
#     if(class(try_result) == "try-error") {
#       attr(res, "msg" <- "cl is not a cluster")
#       return(res)
#     }
#   }
  
  ncl <- length(cl)
  if(ncl < 2) {
    attr(res, "msg") <- "Cluster size is too small (or NULL)"
  } else if (ncl > nboot) {
    attr(res, "msg") <- "nboot is too small for cluster size"
  } else {
    res <- TRUE
  }
  
  return(res)
}
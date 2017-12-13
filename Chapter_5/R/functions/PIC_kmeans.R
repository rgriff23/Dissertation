# This function accepts GPA aligned 3D landmark coordinates [p x 3 x n] and executes the following operations:
# 1: compute phylogenetically independent contrasts (PIC) for each landmark coordinate [p x 3 x n-1]
# 2: use PICs to compute a phylogenetically independent landmark covariance matrix (LCM) [p x p]
# 3: perform k-means clustering for k = 1:p-1
# 4: perform Ward's clustering for k = 1:p-1
# 5: return a list with the following elements
  ## (1) landmark PICs
  ## (2) phylogenetically independent LCM
  ## (3) k-means clustering results
  ## (4) Ward's clustering results 
# The 'ape' package must be loaded 

# function used internally to create list of landmarks in each cluster
#clusters <- function(landmarks_names, cluster_ids) {
#  l <- list()
#  for (i in 1:max(cluster_ids)) {l[[i]] <- landmarks_names[cluster_ids==i]}
#  return(l)
#}

# function to compute PIC-corrected landmark covariance matrix and distance matrix
PIC_cov_dist <- function (coords, phy, allometry=NULL, goswami=FALSE) {
  
  # compute PICs for landmark coordinates - a [p x 3 x n-1] array 
  coords2d <- two.d.array(coords)
  pics <- apply(coords2d, 2, function(x) pic(x, phy=phy))
  
  # control for allometry if allometry != NULL
  if(!is.null(allometry)) {
    allometry_pics <- pic(allometry, phy=phy)
    pics <- apply(pics, 2, function(x) lm(x~allometry_pics)$residuals)
  }
  
  # compute landmark covariance matrix from pics - a [p x p] matrix
  lcm <- matrix(0,dim(coords)[1],dim(coords)[1])
  for (i in 1:dim(coords)[1]) {
    for (j in 1:i) {
      x <- pics[,(3*i-2):(3*i)]
      y <- pics[,(3*j-2):(3*j)]
      if (goswami == FALSE) {
        xy.cov <- cov(cbind(x,y))[1:3,4:6] # xy.cov[1,1] = sum((x[,1] - mean(x[,1]))*(y[,1] - mean(y[,1])))/(n-1)
        xy.cov <- xy.cov*(dim(pics)[1]-1) # get rid of (n-1) 
        covv <- sum(xy.cov^2) # equivalent to sum(diag(xy.cov%*%t(xy.cov)))
        lcm[i,j] <- lcm[j,i] <- covv
      } else lcm[i,j] <- lcm[j,i] <- sum(x*y) # computes goswami-style cov
    }
  }
  dimnames(lcm) <- list(dimnames(coords)[[1]], dimnames(coords)[[1]])
  
  res <- hcut(dist(lcm), k = 6, stand = TRUE)
  fviz_dend(res, rect = TRUE)
  
  # make distance matrix from lcm
  d <- matrix(0,nrow(lcm),nrow(lcm))
  for (i in 2:nrow(lcm)) {
    for (j in 1:(i-1)) {
      euclid <- sqrt(lcm[i,i] + lcm[j,j] - lcm[i,j]*2)
      d[i,j] <- d[j,i] <- euclid
    }
  }
  dimnames(d) <- dimnames(lcm)
  d <- as.dist(d)
  
  # create multidimensional scaling data from the distance matrix
  mds_dat <- suppressWarnings(cmdscale(d, k=(nrow(lcm)-1)))
  mds_dat <- mds_dat[,1:ncol(mds_dat)]
  

  
  # use multidimensional scaling to create data matrix from distance matrix
#  dist <- matrix(0,dim(coords)[1],dim(coords)[1])
#  for (i in 1:(dim(coords)[1]-1)) {
#    for (j in (i+1):dim(coords)[1]) {
#      d2 <- lcm[i,i] + lcm[j,j] - 2*lcm[i,j]
#      dist[i,j] <- dist[j,i] <- sqrt(d2)
#    }
#  }
#  dist <- as.dist(dist)
#  dat <- suppressWarnings(cmdscale(dist,k=dim(coords)[1]-1))
  
  # perform k-means clustering for k = 1:p-1 and save as list of p-1 vectors of length p
#  wclust <- list()
#  for (i in 1:(dim(coords)[1]-1)) {
#    k <- kmeans(dat, i)
#    w <- hcut(dist, k=i)
#    kclust[[i]] <- clusters(dimnames(coords)[[1]], k$cluster)
#    wclust[[i]] <- clusters(dimnames(coords)[[1]], w$cluster)
#  }
  
  # return landmark covariance matrix
  return(list(pics=pics, lcm=lcm, dmat=d, mds=mds_dat))
}

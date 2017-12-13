# This function accepts GPA aligned 3D landmark coordinates [p x 3 x n] 
# and returns a list with the following elements:
  ## (1) phylogenetically independent contrasts (pics) of landmark coordinates
  ## (2) landmark covariance matrix (lcm) derived from pics
  ## (3) Euclidean distance matrix (dmat) derived from lcm
  ## (4) data matrix derived from multidimensional scaling of dmat
# The 'ape' package must be loaded 
require("ape")

# function to compute PIC-corrected landmark covariance matrix and distance matrix
PIC_mats <- function (coords, phy, allometry=NULL, goswami=FALSE) {
  
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
  
  # return landmark covariance matrix
  return(list(pics=pics, lcm=lcm, dmat=d, mds=mds_dat))
}

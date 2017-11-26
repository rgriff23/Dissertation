# This function accepts GPA aligned 3D landmark coordinates [p x 3 x n] and executes the following operations:
# 1: compute phylogenetically independent contrasts (PIC) for each landmark coordinate [p x 3 x n-1]
# 2: use PICs to compute a phylogenetically independent landmark correlation matrix (LCM) [p x p]
# 3: perform k-means clustering for k = 1:p-1
# 4: compute sum of squared errors (SSE) for each value of k 
# 5: return a list with the following elements
  ## (1) landmark PICs
  ## (2) phylogenetically independent LCM
  ## (3) k-means clustering results
  ## (4) SSE for each value of k 
# The 'ape' package must be loaded 
PIC_kmeans <- function (coords, phy) {
  
  # compute PICs for each landmark coordinate and save as [p x 3 x n-1] array 
  coords2d <- two.d.array(coords)
  pics <- apply(coords2d, 2, function(x) pic(x, phy=phy))
  
  # compute phylogenetically independent LCM and save as a [p x p] matrix
  lcm <- matrix(1,dim(coords)[1],dim(coords)[1])
  for (i in 2:dim(coords)[1]) {
    for (j in 1:(i-1)) {
      x <- coords2d[,(3*i-3):(3*i)]
      y <- coords2d[,(3*j-3):(3*j)]
      XY.vcv <- cov(cbind(x,y))
      S12 <- XY.vcv[1:dim(x)[2],(dim(x)[2]+1):(dim(x)[2]+dim(y)[2])]
      S11 <- XY.vcv[1:dim(x)[2],1:dim(x)[2]]
      S22 <- XY.vcv[(dim(x)[2]+1):(dim(x)[2]+dim(y)[2]),(dim(x)[2]+1):(dim(x)[2]+dim(y)[2])]
      rv <- sum(diag(S12%*%t(S12)))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22))) 
      lcm[i,j] <- lcm[j,i] <- rv
    }
  }
  
  # perform k-means clustering for k = 1:p-1 and save as a list of p-1 vectors of length p
  # save SSE for each value of k as a vector of length p-1
  dist <- 1 - lcm
  clusters <- list()
  sse <- c()
  for (i in 1:(dim(coords)[1]-1)) {
    k <- kmeans(dist, i)
    clusters[[i]] <- k$cluster
    sse[i] <- k$tot.withinss
  }
  
  # return list
  results <- list(PICs=pics, LCM=lcm, clust=clusters, SSE=sse)
  return(results)
}
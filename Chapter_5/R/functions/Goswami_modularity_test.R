# This function performs Goswami's modularity test on the results of Ward's cluster analysis
# lcm = a landmark covariance matrix, k = the desired number of clusters
Goswami_modularity_test <- function (lcm, k=6) {
  
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
  
  # Ward's clustering
  res <- hcut(d, k = k)
  clust <- res$cluster
  sizes <- table(res$cluster)
  
  # create correlation matrix (RV coefficient) from lcm
  cm <- matrix(1, nrow(lcm), nrow(lcm))
  for (i in 2:nrow(cm)) {
    for (j in 1:(i-1)) {
      cm[i,j] <- cm[j,i] <- lcm[i,j]/sqrt(lcm[i,i]*lcm[j,j])
    }
  }
  
  # divide correlation matrix into within vs between
  within <- c()
  between <- c()
  for (i in 2:nrow(cm)) {
    for (j in 1:(i-1)) {
      if (clust[i] == clust[j]) {
        within <- c(within, cm[i,j])
      } else between <- c(between, cm[i,j])
    }
  }

  # t-test
  p <- t.test(within, between, alternative="greater")$p.value
  p_out <- noquote(paste("p =", round(p, 4)))
  if (p < 0.05) p_out <- noquote(paste("p =", round(p, 4), "*"))
  if (p < 0.01) p_out <- noquote(paste("p =", round(p, 4), "**"))
  if (p < 0.001) p_out <- noquote(paste("p < 0.001 ***"))
  
  # return
  return(p_out)
}
# Goswami modularity test
Goswami_modularity <- function (lcm, k=6) {
  
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
  
  # divide matrix into within vs between
  within <- c()
  between <- c()
  for (i in 2:nrow(lcm)) {
    for (j in 1:(i-1)) {
      if (clust[i] == clust[j]) {
        within <- c(within, lcm[i,j])
      } else between <- c(between, lcm[i,j])
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
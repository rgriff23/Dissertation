# J. Claude 2008	
# modified from geomorph code

# function to compute centroid size
csize <- function(A) {
  p <- dim(A)[1]
  size <- sqrt(sum(apply(A,2,var))*(p-1))
  list("centroid_size"=size,"scaled"=A/size)
}

# function to center matrix
trans <- function (A) {scale(A, scale=FALSE)}

pPsup <- function (M1,M2) {
  k <- ncol(M1)
  Z1 <- trans(csize(M1)[[2]]) # scale and center
  Z2 <- trans(csize(M2)[[2]]) # scale and center
  sv <- svd(t(Z2)%*%Z1)
  U <- sv$v
  V <- sv$u
  Delt<-sv$d
  sig <- sign(det(t(Z2)%*%Z1))
  Delt[k] <- sig*abs(Delt[k])
  V[,k] <- sig*V[,k]
  Gam <- U%*%t(V)
  # beta<-sum(Delt)   #commented out: retain size-scaling (DCA)
  #  list(Mp1=beta*Z1%*%Gam,Mp2=Z2,rotation=Gam,scale=beta,
  #       df=sqrt(1-beta^2))}
  list(Mp1=Z1%*%Gam, Mp2=Z2, rotation=Gam)} 


pgpa <- function (A) {
  
  # get array dimensions
  p <- dim(A)[1]  # landmarks
  k <- dim(A)[2]  # dimensions
  n <- dim(A)[3]  # configurations
  
  temp1 <- temp2 <- array(NA,dim=c(p,k,n))
  Siz <- numeric(n) # vector of length n to save centroid sizes
  for (i in 1:n) {  # loop through each configuration
    Acs <- csize(A[,,i])  # compute centroid size of configuration
    Siz[i] <- Acs[[1]]  # save centroid size
    temp1[,,i] <- trans(Acs[[2]])}  # center the scaled configuration, add to temp1
    Qm1 <- dist(t(matrix(temp1,k*p,n))) 
    Q <- sum(Qm1)
    iter <- 0
    while (abs(Q) > 0.0001) {# do until change in Q is really tiny
      for (i in 1:n) {
        M <- mshape(temp1[,,-i]) # compute mean shape of all other configurations
        temp2[,,i] <- pPsup(temp1[,,i],M)[[1]]  # 
        }
      Qm2 <- dist(t(matrix(temp2,k*p,n)))
      Q <- sum(Qm1) - sum(Qm2) # change in Q
      Qm1 <- Qm2
      iter <- iter + 1
      temp1 <- temp2
    }
list("rotated"=temp2,
     "it.number"=iter,
     "Q"=Q,
     "intereucl.dist"=Qm2,
     "mshape"= csize(mshape(temp2))[[2]],
     "cent.size"=Siz)
}
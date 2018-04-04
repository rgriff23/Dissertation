# takes 3D coordinate array and returns standard deviation for each coordinate
sd.coords <- function(coords) {
	sds <- matrix(NA,0,3)
	for (l in 1:dim(coords)[1]) {
		sds <- rbind(sds, c(sd(coords[l,1,]),sd(coords[l,2,]),sd(coords[l,3,])))
	}
	rownames(sds) <- dimnames(coords)[[1]]
	colnames(sds) <- c("x","y","z")
	return(sds)
}
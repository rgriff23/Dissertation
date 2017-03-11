# function to read in pp files exported from Meshlab's PickPoints module
# dim assumes 3D
# round=TRUE rounds landmarks to 2 decimal places
read.pp <- function (file, dim=3, round=TRUE) {
	file <- readLines(file)
	lines <- file[grep("point", file)]
	x <- strsplit(lines, "x=\"")
	y <- strsplit(lines, "y=\"")
	z <- strsplit(lines, "z=\"")
	name <- strsplit(lines, "name=\"")
	mat <- matrix(0, length(x), dim)
	r <- c()
	for (i in 1:length(lines)) {
		mat[i,1] <- round(as.numeric(strsplit(x[[i]][2], "\"")[[1]][1]),2)
		mat[i,2] <- round(as.numeric(strsplit(y[[i]][2], "\"")[[1]][1]),2)
		mat[i,3] <- round(as.numeric(strsplit(z[[i]][2], "\"")[[1]][1]),2)
		r <- c(r, strsplit(name[[i]][2], "\"")[[1]][1])
	}
	rownames(mat) <- r 
	colnames(mat) <- c("x","y","z")
	return(mat)
}
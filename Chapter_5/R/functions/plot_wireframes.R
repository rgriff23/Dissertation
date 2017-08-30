# W is a w x 2 matrix with each row representing the start and end points of a wire

# plot wireframe from an n x 3 matrix A and a W matrix
plot.coords <- function(A, W, points.col="black", points.cex=1, lines.col="black", lines.wd=2, bg.col="white", main="", main.line=2, main.cex=2.5, legend=NULL, legend.pos="topright", legend.title="", legend.col=NULL, legend.cex=1.2, legend.lwd=2, legend.bty="n", params=NULL, add=FALSE) {
  if (!is.null(params)) {par3d(params)}
  points.col <- rep(points.col, length.out=nrow(A))
  points.cex <- rep(points.cex, length.out=nrow(A))
  lines.col <- rep(lines.col, length.out=nrow(W))
  lines.wd <- rep(lines.wd, length.out=nrow(W))
  rgl.bg(sphere=TRUE, color=bg.col, lit=FALSE, back="fill")
  plot3d(A, type="s", col=points.col, xlab="", ylab="", zlab="", size=points.cex, aspect=FALSE, box=FALSE, axes=FALSE, add=add)
    if (!is.null(main) | !is.null(legend)) {
      if (!is.null(legend) & is.null(legend.col)) stop("If legend is supplied, must supply legend colors")
      bgplot3d({plot.new()
    if (!is.null(main)) title(main=main, line=main.line, cex.main=main.cex)
    if (!is.null(legend)) legend(legend.pos, title=legend.title, legend=legend, col=legend.col, lwd=legend.lwd, cex=legend.cex, bty=legend.bty)})}
  for (i in 1:nrow(W)) {
    segments3d(rbind(A[W[i,1],], A[W[i,2],]), lwd=lines.wd[i], col=lines.col[i])
  }
}

# plot wireframe from a procD object and a W matrix
# 'value' is the desired value of the last covariate in the model
# the value of all other covariates is assumed to be their average 
plot.procD <- function(procd, W, value=NULL, points.col="black", points.cex=1, lines.col="black", lines.wd=2, bg.col="white", params=NULL, add=FALSE) {
  coeff <- procd$pgls.coefficients
  means <- colMeans(procd$X)
  means[nrow(coeff)] <- value
  coeff <- means*coeff
  A <- matrix(colSums(coeff), ncol(coeff)/3, 3, byrow=TRUE)
  plotCoords(A, W, points.col=points.col, points.cex=points.cex, lines.col=lines.col, lines.wd=lines.wd, bg.col=bg.col, params=params, add=add)
}
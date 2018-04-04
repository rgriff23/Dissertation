# plot wireframe from landmark coordinate matrix
# A is an n x 3 matrix of landmark coordinates
# W is a w x 2 matrix with each row representing the start and end points of a wire
plot.coords <- function(A, W, points.col="black", points.cex=1, lines.col="black", lines.wd=2, bg.col=NULL, 
                        main=NULL, main.line=2, main.cex=2, legend=NULL, legend.pos="topright", legend.title="", 
                        legend.col=NULL, legend.cex=1.2, legend.lwd=2, legend.bty="n", params=NULL, add=FALSE) {
  if (!is.null(params)) {par3d(params)}
  points.col <- rep(points.col, length.out=nrow(A))
  points.cex <- rep(points.cex, length.out=nrow(A))
  lines.col <- rep(lines.col, length.out=nrow(W))
  lines.wd <- rep(lines.wd, length.out=nrow(W))
  if (!is.null(bg.col)) rgl.bg(sphere=TRUE, color=bg.col, lit=FALSE, back="fill")
  plot3d(A, type="s", col=points.col, xlab="", ylab="", zlab="", size=points.cex, aspect=FALSE, box=FALSE, axes=FALSE, add=add)
    if (!is.null(main) | !is.null(legend)) {
      if (!is.null(legend) & is.null(legend.col)) stop("must supply legend colors")
      bgplot3d({plot.new()
    if (!is.null(main)) title(main=main, line=main.line, cex.main=main.cex)
    if (!is.null(legend)) legend(legend.pos, title=legend.title, legend=legend, col=legend.col, lwd=legend.lwd, cex=legend.cex, bty=legend.bty)})}
  for (i in 1:nrow(W)) {
    segments3d(rbind(A[W[i,1],], A[W[i,2],]), lwd=lines.wd[i], col=lines.col[i])
  }
}

# plot wireframe from a procD object and a W matrix
# 'value' is the desired value of the last covariate in the model
# the value of all other covariates is assumed to be their average if means=NULL
# if means != NULL, then means should be a vector of length equal to the number of predictor variables - 1
plot.procD <- function(procd, W, value=NULL, means=NULL, points.col="black", points.cex=1, lines.col="black", lines.wd=2, 
                       bg.col=NULL, main=NULL, main.line=2, main.cex=2, legend=NULL, legend.pos="topright", legend.title="", 
                       legend.col=NULL, legend.cex=1.2, legend.lwd=2, legend.bty="n", params=NULL, add=FALSE) {
  coeff <- procd$pgls.coefficients
  if (is.null(means)) {
    means <- colMeans(procd$X)
    means[nrow(coeff)] <- value
  } else means <- c(1, means, value)
  coeff <- means*coeff
  A <- matrix(colSums(coeff), ncol(coeff)/3, 3, byrow=TRUE)
  plot.coords(A, W, points.col=points.col, points.cex=points.cex, lines.col=lines.col, lines.wd=lines.wd, bg.col=bg.col, 
              main=main, main.line=main.line, main.cex=main.cex, legend=legend, legend.pos=legend.pos, legend.title=legend.title, 
              legend.col=legend.col, legend.cex=legend.cex, legend.bty=legend.bty, params=params, add=add)
}
# A = array containing matrices of of landmarks
# links = matrix with start and end points of lines
# bg_col = background color, link_col = link colors, point_col = point colors
# required packages: rgl
plotWireframe = function (A, links, bg_col="black", link_col="green", point_col="white", params=NULL) 
{
  open3d()
  if (!is.null(params)) {par3d(params)}
  rgl.bg(sphere=TRUE, color=bg_col, lit=FALSE, back="fill")
  A3d = matrix(A, dim(A)[[1]], dim(A)[[2]], dimnames=dimnames(A)[1:2])
  plot3d(A3d, type = "s", col = point_col, xlab = "", ylab = "", zlab = "", size = 1, aspect = FALSE, box=FALSE, axes=FALSE)
  points3d(A, color = point_col, size = 4)
  link_col <- rep(link_col, length.out=nrow(links))
  for (i in 1:nrow(links)) {
    segments3d(rbind(A[links[i,1],,], A[links[i,2],,]), lwd = 2, col=link_col[i])
  }
}


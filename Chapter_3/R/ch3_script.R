################
# PREPARATIONS #
################

# load packages
library("ggplot2")
library("geomorph")
library("scatterplot3d")
library("abind")
library("ggbiplot")
library("gridExtra")

# read landmark data
path <- "Chapter 3 - Jaw closing/Landmarks/"
files <- paste(path, list.files(path=path, pattern=".txt"), sep="")
landmarks <- NULL
for (i in 1:length(files)) {landmarks <- rbind(landmarks, as.matrix(read.table(files[i])))}
landmarks <- arrayspecs(landmarks, 6, 3)

# extract GPA coordinates and centroid sizes
landmarks <- gpagen(landmarks)$coords
sizes <- gpagen(landmarks)$Csize

# subset to mouse and tamarin
mouse <- landmarks[,,1:30]
tamarin <- landmarks[,,31:60]

##################################################################################
# ANALYSIS 1: Compare variances of ground truth, 4-point, and 8-point alignments #
##################################################################################

# GPA landmarks for each specimen
mou_abc <- gpagen(mouse)$coords
tam_abc <- gpagen(tamarin)$coords

# define groups
group30 <- factor(c(rep("ground truth",10),rep("4-point",10),rep("8-point",10)))

# tests of homogenous variance (4-point, 8-point, ground truth)
morphol.disparity(mou_abc ~ group30)
morphol.disparity(tam_abc ~ group30)

############################################################################
# ANALYSIS 2: Compare distances of 4 and 8-point alignment to target shape #
############################################################################

# GPA for landmarks from n-point alignments
tam_bc <- gpagen(tamarin[,,11:30])$coords
mou_bc <- gpagen(mouse[,,11:30])$coords

# compute target shapes
tam_target <- gpagen(tamarin[,,1:10])$consensus
mou_target <- gpagen(mouse[,,1:10])$consensus

# compute Procrustes distances to target for tamarins (20)
tam.dist <- mou.dist <- c()
for (i in 1:20) {
  gpa1 <- gpagen(abind(tam_bc[,,i], tam_target, along=3))$coords
  gpa2 <- gpagen(abind(mou_bc[,,i], mou_target, along=3))$coords
  tam.dist <- c(tam.dist, sum((gpa1[,,1] - gpa1[,,2])^2))
  mou.dist <- c(mou.dist, sum((gpa2[,,1] - gpa2[,,2])^2))
}

# anova (distance ~ method, for each specimen)
group20 <- factor(c(rep("4-point",10),rep("8-point",10)))
anova(lm(tam.dist ~ group20))
anova(lm(mou.dist ~ group20))

# boxplots of distance from target 
boxdat <- data.frame(
  dist=c(mou.dist, tam.dist),
  Method=factor(c(rep(rep(c("4-point", "8-point"), each=10),2))),
  specimen=factor(c(rep(c("Microcebus murinus", "Saguinus fuscicollis"), each=20))))
fig3 <- ggplot(boxdat, aes(specimen, dist, fill=Method)) +
  geom_boxplot(alpha=0.85) + 
  theme(axis.text.x=element_text(face=4, size = 10),
        axis.text.y=element_text(size=11),
        axis.title.y=element_text(margin=margin(0,20,0,0), size = 13, face = 2),
        legend.title=element_text(size=13),
        legend.text=element_text(size=11)) +
  ylab("Procrustes distance from ground truth") + xlab("")  +
  scale_fill_manual(values=c("red4","hotpink2"))
fig3

#############################################
# ANALYSIS 3: Principal components analysis #
#############################################

# pca
tam_2d <- two.d.array(tam_abc)
mou_2d <- two.d.array(mou_abc)
all_2d <- two.d.array(landmarks)
colnames(tam_2d) <- colnames(mou_2d) <- colnames(all_2d) <- as.character(1:18)
tam.pca <- prcomp(tam_2d)
mou.pca <- prcomp(mou_2d)
all.pca <- prcomp(all_2d)

# create pca min and max shapes for tamarin
min1 <- max1 <- rep(0,18)
min1[1] <- min(tam.pca$x[,1])
max1[1] <- max(tam.pca$x[,1])
tam.min <- arrayspecs(as.matrix(min1 %*% (t(tam.pca$rotation))),6,3)[,,1] + mshape(tam_abc)
tam.max <- arrayspecs(as.matrix(max1 %*% (t(tam.pca$rotation))),6,3)[,,1] + mshape(tam_abc)

# create pca min and max shapes for mouse lemur
min1[1] <- min(mou.pca$x[,1])
max1[1] <- max(mou.pca$x[,1])
mou.min <- arrayspecs(as.matrix(min1 %*% (t(mou.pca$rotation))),6,3)[,,1] + mshape(mou_abc)
mou.max <- arrayspecs(as.matrix(max1 %*% (t(mou.pca$rotation))),6,3)[,,1] + mshape(mou_abc)

# barplots showing variance across pcas
bardat <- c(summary(mou.pca)$importance[2,], summary(tam.pca)$importance[2,])
bardat <- data.frame(bars=bardat,
                     order=factor(names(bardat), levels=names(bardat)[1:18]),
                     specimen=factor(c(rep(c("Microcebus murinus", "Saguinus fuscicollis"), each=18))))
ggplot(bardat, aes(order, bars)) +
  geom_bar(stat="identity") +
  facet_grid( ~ specimen) +
  ylab("Proportion of variance explained") + xlab("Principle components") +
  theme(axis.text.x=element_text(angle=90),
        axis.title.x=element_text(face=2, size=13),
        axis.title.y=element_text(margin=margin(0,20,0,0), size=13, face=2),
        strip.text = element_text(face = 3, size=13))

# biplots for PCA of each specimen separately
groups3 <- rep(c("Ground truth", "4-point", "8-point"), each=10)
ggbi1 <- ggbiplot(mou.pca, groups=groups3, ellipse=TRUE, var.axes=FALSE) +
  labs(color="Method") +
  scale_color_manual(values=c("red4","hotpink2","gold3")) +
  ggtitle("Microcebus murinus") +
  theme(axis.title = element_text(size=13),
        plot.title=element_text(face="italic", hjust=0.5))
ggbi2 <- ggbiplot(tam.pca, groups=groups3, ellipse=TRUE, var.axes=FALSE) +
  labs(color="Method") +
  scale_color_manual(values=c("blue4","dodgerblue","purple")) +
  ggtitle("Saguinus fuscicollis") +
  theme(axis.title = element_text(size=13),
        plot.title=element_text(face="italic", hjust=0.5))
grid.arrange(ggbi1, ggbi2, ncol=2)

# visualize shape changes along PC1
layout3d(matrix(1:2, 1, 2),sharedMouse = FALSE)
plot3d(c(mou.min[,1], mou.max[,1]), c(mou.min[,2], mou.max[,2]), c(mou.min[,3], mou.max[,3]), 
       xlab="", ylab="", zlab="", col=rep(c("gold3","hotpink2"), each=6), size=5, box=F, axes=F)
bgplot3d({plot.new()
  title(main=expression(italic("Microcebus murinus")), line=2, cex.main=2)
  legend("bottomright", legend=c("Physical", "n-point"), title="Alignment type", col=c("gold3", "hotpink2"), pch=16, lwd=2, cex=1.5, bty="n")})
segments3d(mou.min, col="gold3")
lines3d(mou.min[c(1,3,5,1,2,4,6,2),], col="gold3", lwd=3)
triangles3d(mou.min[c(1,3,5,2,4,6),],col="gold3")
lines3d(mou.max[c(1,3,5,1,2,4,3,4,6,2,6,5),], col="hotpink2", lwd=3)
plot3d(c(tam.min[,1], tam.max[,1]), c(tam.min[,2], tam.max[,2]), c(tam.min[,3], tam.max[,3]), 
       xlab="", ylab="", zlab="", col=rep(c("gold3","hotpink2"), each=6), size=5, box=F, axes=F)
bgplot3d({plot.new()
  title(main=expression(italic("Saguinus fuscicollis")), line=2, cex.main=2)})
segments3d(tam.min, col="gold3")
lines3d(tam.min[c(1,3,5,1,2,4,6,2),], col="gold3", lwd=3)
triangles3d(tam.min[c(1,3,5,2,4,6),],col="gold3")
lines3d(tam.max[c(1,3,5,1,2,4,3,4,6,2,6,5),], col="hotpink2", lwd=3)
#legend3d("topleft", legend=c("Target shape", "8-point alignment"), col=c("blue", "red4"), pch=16)

# biplot for PCA of all specimens
all.group <- factor(c(rep("Microcebus, ground truth",10),
                      rep("Microcebus, 4-point",10),
                      rep("Microcebus, 8-point",10),
                      rep("Saguinus, ground truth",10),
                      rep("Saguinus, 4-point",10),
                      rep("Saguinus, 8-point",10)))
ggbiplot(all.pca, groups= all.group, ellipse=TRUE, var.axes=FALSE, alpha=0.5, scale=0) +
  labs(color="Method") +
  scale_color_manual(values=c("red4","hotpink2","gold3","blue4","dodgerblue","purple")) +
  theme(axis.title = element_text(size=12),
        legend.title = element_text(size=13),
        legend.text = element_text(size=10)) 

##############################################
# Analysis 4: Absolute vertical displacement #
##############################################

# absolute (mm) displacement of cranium and mandible under 8-point alignment relative to ground truth
mou_molar_dist8 <- mean(dist[dist$Pair %in% c("1-2", "5-6") & dist$Specimen == "Microcebus murinus",][1:20,"Distance"]) -
  mean(dist[dist$Pair %in% c("1-2", "5-6") & dist$Specimen == "Microcebus murinus",][41:60,"Distance"])
tam_molar_dist8 <- mean(dist[dist$Pair %in% c("1-2", "5-6") & dist$Specimen == "Saguinus fuscicollis",][1:20,"Distance"]) -
  mean(dist[dist$Pair %in% c("1-2", "5-6") & dist$Specimen == "Saguinus fuscicollis",][41:60,"Distance"])
mou_incisor_dist8 <- mean(dist[dist$Pair %in% c("3-4") & dist$Specimen == "Microcebus murinus",][1:10,"Distance"]) -
  mean(dist[dist$Pair %in% c("3-4") & dist$Specimen == "Microcebus murinus",][21:30,"Distance"])
tam_incisor_dist8 <- mean(dist[dist$Pair %in% c("3-4") & dist$Specimen == "Saguinus fuscicollis",][1:10,"Distance"]) -
  mean(dist[dist$Pair %in% c("3-4") & dist$Specimen == "Saguinus fuscicollis",][21:30,"Distance"])
mou_molar_dist8
tam_molar_dist8
mou_incisor_dist8
tam_incisor_dist8
mou_molar_dist8/31
tam_molar_dist8/45

# absolute (mm) displacement of cranium and mandible under 4-point alignment relative to ground truth
mou_molar_dist4 <- mean(dist[dist$Pair %in% c("1-2", "5-6") & dist$Specimen == "Microcebus murinus",][1:20,"Distance"]) -
  mean(dist[dist$Pair %in% c("1-2", "5-6") & dist$Specimen == "Microcebus murinus",][21:40,"Distance"])
tam_molar_dist4 <- mean(dist[dist$Pair %in% c("1-2", "5-6") & dist$Specimen == "Saguinus fuscicollis",][1:20,"Distance"]) -
  mean(dist[dist$Pair %in% c("1-2", "5-6") & dist$Specimen == "Saguinus fuscicollis",][21:40,"Distance"])
mou_incisor_dist4 <- mean(dist[dist$Pair %in% c("3-4") & dist$Specimen == "Microcebus murinus",][1:10,"Distance"]) -
  mean(dist[dist$Pair %in% c("3-4") & dist$Specimen == "Microcebus murinus",][11:20,"Distance"])
tam_incisor_dist4 <- mean(dist[dist$Pair %in% c("3-4") & dist$Specimen == "Saguinus fuscicollis",][1:10,"Distance"]) -
  mean(dist[dist$Pair %in% c("3-4") & dist$Specimen == "Saguinus fuscicollis",][11:20,"Distance"])
mou_molar_dist4
tam_molar_dist4
mou_incisor_dist4
tam_incisor_dist4
mou_molar_dist4/31
tam_molar_dist4/45

# distances between pairs of points in mm (1-2, 3-4, and 5-6 for each of the two specimens)
dist <- data.frame(
  Specimen=factor(c(rep(c("Microcebus murinus", "Saguinus fuscicollis"), each=90))),
  Method=factor(rep(rep(c("Ground truth","4-point","8-point"), each=30),2)),
  Pair=factor(rep(c("1-2","3-4","5-6"), 60))
)
Distance <- c()
all_2d_mm <- all_2d*sizes
for (i in 1:60) {
  ABC <- c()
  ABC[1] <- sqrt(sum((all_2d_mm[i,1:3] - all_2d_mm[i,4:6])^2))
  ABC[2] <- sqrt(sum((all_2d_mm[i,7:9] - all_2d_mm[i,10:12])^2))
  ABC[3] <- sqrt(sum((all_2d_mm[i,13:15] - all_2d_mm[i,16:18])^2))
  Distance <- c(Distance, ABC)
}
dist$Distance <- Distance
ggplot(dist, aes(Pair, Distance, fill=Method)) +
  geom_boxplot() +
  facet_grid( ~ Specimen) +
  xlab("Pairs of landmarks") + ylab("Vertical distance between landmarks (mm)") +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(margin=margin(20,0,0,0), face=2, size=13),
        axis.title.y=element_text(margin=margin(0,20,0,0), size=13, face=2),
        strip.text = element_text(face = 3, size=13)) +
  scale_fill_manual(values=c("red4","hotpink2","gold3","blue4","dodgerblue","purple"))

########
# END ##
########
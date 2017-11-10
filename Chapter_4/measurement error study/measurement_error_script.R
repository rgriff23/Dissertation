################
# PREPARATIONS #
################

# load packages
library(geomorph)
library(ggplot2) 
library(plyr) # for 'mapvalues' function
library(abind) # for formatting landmark array

# load custom functions
source('~/Desktop/GitHub/Dissertation/Chapter_4/R/functions/read.pp.R', chdir = TRUE)
source('~/Desktop/GitHub/Dissertation/Chapter_4/R/functions/sd.coords.R', chdir = TRUE)

# read landmarks
path = "~/Desktop/GitHub/Dissertation/Chapter_4/data/"
files <- paste(path, list.files(path=path, pattern=".pp"), sep="")
landmarks <- NULL
for (i in 1:length(files)) {landmarks <- abind(landmarks, read.pp(files[i]), along=3)}

############
# ANALYSIS #
############

# compute standard deviations for each landmark/dimension 
sds.fas <- data.frame(round(sd.coords(landmarks[,,1:16]),3))
sds.nem <- data.frame(round(sd.coords(landmarks[,,17:32]),3))
sds.fas$ave <- round(sqrt(rowSums(sds.fas^2)/3),3)
sds.nem$ave <- round(sqrt(rowSums(sds.nem^2)/3),3)
sds <- cbind(sds.fas, sds.nem)
names(sds) <- c("fas.x","fas.y","fas.z","fas.average","nem.x","nem.y","nem.z","nem.average")
sds

# nested ANOVA
alignment <- factor(rep(rep(1:4, each=4), 2))
specimen <- rep(c("M.fas", "M.nem"), each=16)
landmarks.pca <- gpagen(landmarks)$coords
mod <- procD.lm(landmarks.pca ~ specimen/alignment)

#########
# PLOTS #
#########

# define wireframe, colors, and mean shape
wireframe <- read.csv("~/Desktop/GitHub/Dissertation/Chapter_4/data/wireframe.csv", header=F)
modules <- c("mandible", "face", "braincase", "zygomatic", "basicranium")
module_colors <- c("red", "green", "blue", "purple", "gold")
wireframe$colors <- mapvalues(wireframe[,1], modules, module_colors)
mean.fas <- gpagen(landmarks[,,1:16])$consensus
mean.nem <- gpagen(landmarks[,,17:32])$consensus

# 3D wireframes
layout3d(matrix(1:2, 1, 2), sharedMouse = TRUE)
plot3d(mean.fas, xlab="", ylab="", zlab="", size=1, col="darkgray", type='s', aspect=F, box=F, axes=F)
bgplot3d({plot.new()
  title(main=expression(italic("M. fascicularis")), line=2, cex.main=2.5)
  legend("topright", title="Region", legend=modules, col=module_colors, lwd=2, cex=1.2, bty="n")})
for (i in 1:nrow(wireframe)) {
  segments3d(rbind(mean.fas[wireframe[i,2],], mean.fas[wireframe[i,3],]), lwd = 2, col=wireframe[i,4])
}
plot3d(mean.nem, xlab="", ylab="", zlab="", size=1, col="darkgray", type='s', aspect=F, box=F, axes=F)
bgplot3d({plot.new()
  title(main=expression(italic("M. nemestrina")), line=2, cex.main=2.5)})
for (i in 1:nrow(wireframe)) {
  segments3d(rbind(mean.nem[wireframe[i,2],], mean.nem[wireframe[i,3],]), lwd = 2, col=wireframe[i,4])
}
snapshot3d("~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Landmarks and measurement error/Figures/Figure 4.2 Wireframes.png", fmt="png")

# histograms showing distribution of coordinate standard deviations
hist_dat <- data.frame(sds=c(sds.fas, sds.nem), 
  sp=factor(rep(c("M. fascicularis","M. nemestrina"),each=26*3)),
  v=rep(c(median(sds.fas),median(sds.nem)),each=26*3))
ggplot(hist_dat, aes(sds)) +
  geom_histogram(binwidth=0.05, boundary=0.5, color="white") +
  facet_grid( ~ sp) + ylab("Frequency") + xlab("Standard deviation (mm)") +
  ylim(c(0,25)) +
  theme(axis.title.y=element_text(margin=margin(0,20,0,5), size=13, face=2),
  axis.title.x=element_text(margin=margin(20,0,0,5), size=13, face=1),
  strip.text.x=element_text(face="italic", size=14),
  plot.margin=unit(c(1,1,0.5,0),"cm")) +
  geom_vline(aes(xintercept=v), linetype=2, size=1)

########
# END ##
########





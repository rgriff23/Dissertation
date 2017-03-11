################
# PREPARATIONS #
################

# load packages
library(geomorph)
library(ggplot2)
library(plyr) # for 'mapvalues' function to create wireframe colors
library(abind) # for reading in landmarks

# load custom functions
source('~/Desktop/GitHub/Dissertation/Chapter_4/R/functions/read.pp.R', chdir = TRUE)
source('~/Desktop/GitHub/Dissertation/Chapter_4/R/functions/sd.coords.R', chdir = TRUE)
source('~/Desktop/GitHub/Dissertation/Chapter_4/R/functions/plotWireframe.R', chdir = TRUE)

# read landmarks
path = "~/Desktop/GitHub/Dissertation/Chapter_4/data/"
files <- paste(path, list.files(path=path, pattern=".pp"), sep="")
landmarks <- NULL
for (i in 1:length(files)) {landmarks <- abind(landmarks, read.pp(files[i]), along=3)}

############
# ANALYSIS #
############

# compute standard deviations for each landmark/dimension 
sds.fas <- sd.coords(landmarks[,,1:16])
sds.nem <- sd.coords(landmarks[,,17:32])
sds.fas
sds.nem

# nested ANOVA
alignment <- factor(rep(rep(1:4, each=4), 2))
specimen <- rep(c("M.fas", "M.nem"), each=16)
landmarks.pca <- gpagen(landmarks)$coords
mod <- procD.lm(landmarks.pca ~ specimen/alignment)
mod

#########
# PLOTS #
#########

# distribution of coordinate standard deviations
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

# 3d scatterplot
plotAllSpecimens(landmarks[,,1:16], mean=T, plot.param=list(pt.cex=0.5))

# define wireframe
wireframe <- read.csv("~/Desktop/GitHub/Dissertation/Chapter_4/data/wireframe.csv", header=F)

# define wire colors
wire_colors <- mapvalues(wireframe[,1], c("mandible", "face", "braincase", "zygomatic", "basicranium"), c("red", "green", "blue", "purple", "yellow"))

# import landmark data
test <- array(landmarks[,,1], dim=c(26,3,1))
plotWireframe(test, wireframe[,2:3], link_col=wire_colors, bg_col="white")

########
# END ##
########





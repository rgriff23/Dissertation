################
# PREPARATIONS #
################

# load packages
library(geomorph)
library(plyr) # for 'mapvalues' function to create wireframe colors
library(abind) # for reading in landmarks

# load custom functions
source('~/Desktop/GitHub/Dissertation/Chapter_4/R/functions/read.pp.R', chdir = TRUE)
source('~/Desktop/GitHub/Dissertation/Chapter_4/R/functions/sd.coords.R', chdir = TRUE)
source('~/Desktop/GitHub/Dissertation/Chapter_4/R/functions/plotWireframe.R', chdir = TRUE)

# read landmarks
path = "~/Desktop/GitHub/Dissertation/Chapter_4/data"
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

# 3d scatterplot
plotAllSpecimens(landmarks, mean=T, plot.param=list(pt.cex=0.5))

# define wireframe
wireframe <- read.csv("~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Landmarks and measurement error/Spreadsheets/wireframe.csv", header=F)

# define wire colors
wire_colors <- mapvalues(wireframe[,1], c("mandible", "face", "braincase", "zygomatic", "basicranium"), c("red", "green", "blue", "purple", "yellow"))

# import landmark data
#test <- read.pp("~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Landmarks and measurement error/Landmarks/nem_4_4.pp")
#test <- array(test, dim=c(26,3,1))

# define mean shapes for M.fas and M.nem

# make side-by-side wireframe plots with legend
#plotWireframe(test, wireframe[,2:3], link_col=wire_colors)

########
# END ##
########





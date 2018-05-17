################
# PREPARATIONS #
################

# load packages
library("geomorph")
library("plyr")
library("tidyverse")
library("abind") 
library("caper")

# load custom functions
source("Chapter_4/R/functions/read.pp.R")
source("Chapter_4/R/functions/plot_wireframes.R")

# read & tidy data
data <- read.csv("Chapter_4/data/primate_data.csv") %>%
  mutate(Folivore = ifelse(Diet %in% c("Folivore", "Frug-Fol"), 1, 0), 
         id = paste(Genus, Species, Sex, MuseumID, sep="_"), 
         genus_species = paste(Genus, Species, sep="_")) %>%
  select(1:15, Folivore, id, genus_species) %>%
  arrange(id)

# read phylogeny
tree <- read.nexus("Chapter_2/data/tree.nex")

# read landmarks & add dimnames
path = "Chapter_4/data/"
files <- paste(path, list.files(path=path, pattern=".pp"), sep="")
landmarks <- NULL
for (i in 1:length(files)) {landmarks <- abind(landmarks, read.pp(files[i]), along=3)}
dimnames(landmarks)[[3]] <- unlist(strsplit(unlist(strsplit(files, path))[seq(from=2, to=244, length.out=122)], "_picked_points.pp"))
rm(path, files)

# separate data for males and females
data.f <- filter(data, SexDat %in% c("F","B")) %>% select(-Sex, -SexDat)
data.m <- filter(data, SexDat %in% c("M","B")) %>% select(-Sex, -SexDat)

# separate tree for males and females
tree.f <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% data.f$genus_species)])
tree.m <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% data.m$genus_species)])

# separate landmarks for males and females
landmarks.f <- landmarks[,,dimnames(landmarks)[[3]]%in%data.f$id]
dimnames(landmarks.f)[[3]] <- data.f$genus_species
landmarks.m <- landmarks[,,dimnames(landmarks)[[3]]%in%data.m$id]
dimnames(landmarks.m)[[3]] <- data.m$genus_species

# estimate missing landmarks for males and females
landmarks.f <- estimate.missing(landmarks.f, method="Reg")
dimnames(landmarks.f)[1:2] <- dimnames(landmarks)[1:2]
landmarks.m <- estimate.missing(landmarks.m, method="Reg")
dimnames(landmarks.m)[1:2] <- dimnames(landmarks)[1:2]

# procrustes alignment for males and females
gpa.f <- gpagen(landmarks.f)
gpa.m <- gpagen(landmarks.m)

# procrustes alignment of brains for males and females
brains <- c("Vertex","Opisthocranium","Opisthion","Basion","Floor of sella","Glabella","Euryon","Frontotemporale")
gpa.f2 <- gpagen(landmarks.f[brains,,])
gpa.m2 <- gpagen(landmarks.m[brains,,])

# procrustes alignment of crania for males and females
crania <- c("Rhinion","Prosthion","Canine (maxillary)","M2 (maxillary)","Porion","Zygomatic point","Orbit margin (lateral)","Orbit margin (inferior)","Orbit margin (medial)","Orbit margin (superior)","Glabella","Frontotemporale","Euryon","Vertex","Opisthocranium","Opisthion","Basion","Floor of sella") 
gpa.f3 <- gpagen(landmarks.f[crania,,])
gpa.m3 <- gpagen(landmarks.m[crania,,])

# add brain, crania, and skull centroid size to data
data.f$Cranium_size <- gpa.f3$Csize
data.m$Cranium_size <- gpa.m3$Csize
data.f$Brain_size <- gpa.f2$Csize
data.m$Brain_size <- gpa.m2$Csize
data.f$Skull_size <- gpa.f$Csize
data.m$Skull_size <- gpa.m$Csize

########
# PGLS #
########

# estimate brain ~ skull allometric slope with pgls
compdat.f <- comparative.data(tree.f, data.f, names.col=genus_species)
compdat.m <- comparative.data(tree.m, data.m, names.col=genus_species)
pgls.f <- pgls(log(Brain_size) ~ log(Skull_size), compdat.f)
pgls.m <- pgls(log(Brain_size) ~ log(Skull_size), compdat.m)
pgls.f # coefficient = 0.84
pgls.m # coefficient = 0.75

# estimate brain ~ cranium allometric slope with pgls
compdat.f2 <- comparative.data(tree.f, data.f, names.col=genus_species)
compdat.m2 <- comparative.data(tree.m, data.m, names.col=genus_species)
pgls.f2 <- pgls(log(Brain_size) ~ log(Cranium_size), compdat.f2)
pgls.m2 <- pgls(log(Brain_size) ~ log(Cranium_size), compdat.m2)
pgls.f2 # coefficient = 0.89
pgls.m2 # coefficient = 0.81

###########################################
# non-phylogenetic models for skull shape #
###########################################

# define wireframe plot function for non-phylogenetic procD object
plot.procD.lm <- function(procd, W, value=NULL, means=NULL, points.col="black", points.cex=1, lines.col="black", lines.wd=2, 
                       bg.col=NULL, main=NULL, main.line=2, main.cex=2, legend=NULL, legend.pos="topright", legend.title="", 
                       legend.col=NULL, legend.cex=1.2, legend.lwd=2, legend.bty="n", params=NULL, add=FALSE) {
  coeff <- procd$coefficients
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

# create geomorph dataframe for males and females
gdf.f <- geomorph.data.frame(gpa.f)
gdf.m <- geomorph.data.frame(gpa.m)

# shape vs size
nophy.pgls.f <- procD.lm(coords ~ log(Csize), data=gdf.f) 
nophy.pgls.m <- procD.lm(coords ~ log(Csize), data=gdf.m) 
nophy.pgls.f # p < 0.001 ***
nophy.pgls.m # p < 0.001 ***

# define wireframe
wireframe <- read.csv("./Chapter_4/data/wireframe.csv", header=FALSE)

# wireframe plot
layout3d(matrix(1:4, 2, 2, byrow=TRUE), sharedMouse = TRUE)
plot.procD.lm(nophy.pgls.f, wireframe[,2:3], value=log(min(gdf.f$Csize)))
plot.procD.lm(nophy.pgls.f, wireframe[,2:3], value=log(max(gdf.f$Csize)))
plot.procD.lm(nophy.pgls.m, wireframe[,2:3], value=log(min(gdf.m$Csize)))
plot.procD.lm(nophy.pgls.m, wireframe[,2:3], value=log(max(gdf.m$Csize)))

# same thing, cranium only
gdf.f3 <- geomorph.data.frame(gpa.f3, phy=tree.f)
gdf.m3 <- geomorph.data.frame(gpa.m3, phy=tree.m)

# shape vs size, cranium only, non-phylogenetic
nophy.pgls.f3 <- procD.lm(coords ~ log(Csize), data=gdf.f3) 
nophy.pgls.m3 <- procD.lm(coords ~ log(Csize), data=gdf.m3) 
nophy.pgls.f3 # p < 0.001 ***
nophy.pgls.m3 # p < 0.001 ***

# wireframe
col1 <- c(1,2,3,5,6,7,9,7,9,10,1,7,4,10,5,12,5,5,5,11,14,13,15,16,17,11)
col2 <- c(2,3,4,6,8,8,8,10,10,11,11,12,12,12,12,13,13,16,17,14,15,15,16,17,18,18)
wireframe2 <- cbind(col1,col2)

# wireframe plot (cranium only, non-phylogenetic)
layout3d(matrix(1:4, 2, 2, byrow=TRUE), sharedMouse = TRUE)
plot.procD.lm(nophy.pgls.f3, wireframe2, value=log(min(gdf.f3$Csize)))
plot.procD.lm(nophy.pgls.f3, wireframe2, value=log(max(gdf.f3$Csize)))
plot.procD.lm(nophy.pgls.m3, wireframe2, value=log(min(gdf.m3$Csize)))
plot.procD.lm(nophy.pgls.m3, wireframe2, value=log(max(gdf.m3$Csize)))

########
# END ##
########

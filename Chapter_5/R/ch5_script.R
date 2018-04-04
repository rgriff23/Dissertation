################
# PREPARATIONS #
################

# load packages
library("geomorph")
library("Morpho")
library("tidyverse")
library("plyr")
library("abind") 
library("cluster")
library("factoextra")
library("gridExtra")

# load custom functions
source("Chapter_4/R/functions/read.pp.R")
source("Chapter_4/R/functions/plot_wireframes.R")
source("Chapter_5/R/functions/PIC_mats.R")
source("Chapter_5/R/functions/Goswami_modularity_test.R")

# read & tidy data
data <- read.csv("Chapter_4/data/primate_data.csv") %>%
  mutate(id = paste(Genus, Species, Sex, MuseumID, sep="_"), 
         genus_species = paste(Genus, Species, sep="_")) %>%
  select(Superfamily, Sex, SexDat, DimorphismIndex, id, genus_species) %>%
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

# separate data, tree, and landmarks for males and females
data.f <- filter(data, SexDat %in% c("F","B")) %>% select(-Sex, -SexDat)
data.m <- filter(data, SexDat %in% c("M","B")) %>% select(-Sex, -SexDat)
landmarks.f <- landmarks[,,dimnames(landmarks)[[3]]%in%data.f$id]
dimnames(landmarks.f)[[3]] <- data.f$genus_species
landmarks.m <- landmarks[,,dimnames(landmarks)[[3]]%in%data.m$id]
dimnames(landmarks.m)[[3]] <- data.m$genus_species

# estimate missing landmarks for males and females
landmarks.f <- estimate.missing(landmarks.f, method="Reg")
dimnames(landmarks.f)[1:2] <- dimnames(landmarks)[1:2]
landmarks.m <- estimate.missing(landmarks.m, method="Reg")
dimnames(landmarks.m)[1:2] <- dimnames(landmarks)[1:2]

# list of strepsirrhine/haplorhine males/females
streps.f <- data.f[data.f$Superfamily%in%c("Galagoidea","Lemuroidea"),c("DimorphismIndex", "genus_species")]
streps.m <- data.m[data.m$Superfamily%in%c("Galagoidea","Lemuroidea"),c("DimorphismIndex", "genus_species")]
haps.f <- data.f[data.f$Superfamily%in%c("Cercopithecoidea","Ceboidea","Hominoidea"),c("DimorphismIndex", "genus_species")]
haps.m <- data.m[data.m$Superfamily%in%c("Cercopithecoidea","Ceboidea","Hominoidea"),c("DimorphismIndex", "genus_species")]

# subset tree for male/female strepsirrhines/haplorhines
tree.primates.f <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% data.f$genus_species)])
tree.primates.m <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% data.m$genus_species)])
tree.streps.f <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% streps.f$genus_species])
tree.streps.m <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% streps.m$genus_species])
tree.haps.f <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% haps.f$genus_species])
tree.haps.m <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% haps.m$genus_species])

# procrustes alignment for males and females (all primates, strepsirrhines, haplorhines)
gpa.primates.f <- gpagen(landmarks.f)
gpa.primates.m <- gpagen(landmarks.m)
gpa.streps.f <- gpagen(landmarks.f[,,dimnames(landmarks.f)[[3]]%in%streps.f$genus_species])
gpa.streps.m <- gpagen(landmarks.m[,,dimnames(landmarks.m)[[3]]%in%streps.m$genus_species])
gpa.haps.f <- gpagen(landmarks.f[,,dimnames(landmarks.f)[[3]]%in%haps.f$genus_species])
gpa.haps.m <- gpagen(landmarks.m[,,dimnames(landmarks.m)[[3]]%in%haps.m$genus_species])

# size-preserving Procrustes superimposition
gpa2.primates.f <- procSym(landmarks.f, reflect=FALSE, scale=FALSE, CSinit=FALSE)$rotated
gpa2.primates.m <- procSym(landmarks.m, reflect=FALSE, scale=FALSE, CSinit=FALSE)$rotated
gpa2.streps.f <- procSym(landmarks.f[,,dimnames(landmarks.f)[[3]]%in%streps.f$genus_species], reflect=FALSE, scale=FALSE, CSinit=FALSE)$rotated
gpa2.streps.m <- procSym(landmarks.m[,,dimnames(landmarks.m)[[3]]%in%streps.m$genus_species], reflect=FALSE, scale=FALSE, CSinit=FALSE)$rotated
gpa2.haps.f <- procSym(landmarks.f[,,dimnames(landmarks.f)[[3]]%in%haps.f$genus_species], reflect=FALSE, scale=FALSE, CSinit=FALSE)$rotated
gpa2.haps.m <- procSym(landmarks.m[,,dimnames(landmarks.m)[[3]]%in%haps.m$genus_species], reflect=FALSE, scale=FALSE, CSinit=FALSE)$rotated

# create matrices of control variables (centroid size, dimorphism index) for each dataset
cvar.primates.f <- matrix(c(log(gpa.primates.f$Csize), data.f$DimorphismIndex), ncol=2)
cvar.primates.m <- matrix(c(log(gpa.primates.m$Csize), data.m$DimorphismIndex), ncol=2)
cvar.streps.f <- matrix(c(log(gpa.streps.f$Csize), streps.f$DimorphismIndex), ncol=2)
cvar.streps.m <- matrix(c(log(gpa.streps.m$Csize), streps.m$DimorphismIndex), ncol=2)
cvar.haps.f <- matrix(c(log(gpa.haps.f$Csize), haps.f$DimorphismIndex), ncol=2)
cvar.haps.m <- matrix(c(log(gpa.haps.m$Csize), haps.m$DimorphismIndex), ncol=2)

####################################################################################################
# COMPUTE PIC-BASED LANDMARK VCV MATRICES, EUCLIDEAN DISTANCE MATRICES, & MULTIDIMENSIONAL SCALING #
####################################################################################################

# no control variables 
u.primates.f <- PIC_mats(gpa.primates.f$coords, tree.primates.f)
u.primates.m <- PIC_mats(gpa.primates.m$coords, tree.primates.m)
u.streps.f <- PIC_mats(gpa.streps.f$coords, tree.streps.f)
u.streps.m <- PIC_mats(gpa.streps.m$coords, tree.streps.m)
u.haps.f <- PIC_mats(gpa.haps.f$coords, tree.haps.f)
u.haps.m <- PIC_mats(gpa.haps.m$coords, tree.haps.m)

# include control variables
c.primates.f <- PIC_mats(gpa.primates.f$coords, tree.primates.f, controlVars=cvar.primates.f)
c.primates.m <- PIC_mats(gpa.primates.m$coords, tree.primates.m, controlVars=cvar.primates.m)
c.streps.f <- PIC_mats(gpa.streps.f$coords, tree.streps.f, controlVars=cvar.streps.f)
c.streps.m <- PIC_mats(gpa.streps.m$coords, tree.streps.m, controlVars=cvar.streps.m)
c.haps.f <- PIC_mats(gpa.haps.f$coords, tree.haps.f, controlVars=cvar.haps.f)
c.haps.m <- PIC_mats(gpa.haps.m$coords, tree.haps.m, controlVars=cvar.haps.m)

# Goswami's method
g.primates.f <- PIC_mats(gpa2.primates.f, tree.primates.f, goswami=TRUE, landmark_names=dimnames(landmarks)[[1]])
g.primates.m <- PIC_mats(gpa2.primates.m, tree.primates.m, goswami=TRUE, landmark_names=dimnames(landmarks)[[1]])
g.streps.f <- PIC_mats(gpa2.streps.f, tree.streps.f, goswami=TRUE, landmark_names=dimnames(landmarks)[[1]])
g.streps.m <- PIC_mats(gpa2.streps.m, tree.streps.m, goswami=TRUE, landmark_names=dimnames(landmarks)[[1]])
g.haps.f <- PIC_mats(gpa2.haps.f, tree.haps.f, goswami=TRUE, landmark_names=dimnames(landmarks)[[1]])
g.haps.m <- PIC_mats(gpa2.haps.m, tree.haps.m, goswami=TRUE, landmark_names=dimnames(landmarks)[[1]])

######################################
# CLUSTER ANALYSIS AND VISUALIZATION #
######################################

# Ward's clustering 
k <- 6
# no control variables
res.u.primates.f <- hcut(u.primates.f$dmat, k=k)
res.u.primates.m <- hcut(u.primates.m$dmat, k=k)
res.u.streps.f <- hcut(u.streps.f$dmat, k=k)
res.u.streps.m <- hcut(u.streps.m$dmat, k=k)
res.u.haps.f <- hcut(u.haps.f$dmat, k=k)
res.u.haps.m <- hcut(u.haps.m$dmat, k=k)
# include control variables
res.c.primates.f <- hcut(c.primates.f$dmat, k=k)
res.c.primates.m <- hcut(c.primates.m$dmat, k=k)
res.c.streps.f <- hcut(c.streps.f$dmat, k=k)
res.c.streps.m <- hcut(c.streps.m$dmat, k=k)
res.c.haps.f <- hcut(c.haps.f$dmat, k=k)
res.c.haps.m <- hcut(c.haps.m$dmat, k=k)
# Goswami's method
res.g.primates.f <- hcut(g.primates.f$dmat, k=k)
res.g.primates.m <- hcut(g.primates.m$dmat, k=k)
res.g.streps.f <- hcut(g.streps.f$dmat, k=k)
res.g.streps.m <- hcut(g.streps.m$dmat, k=k)
res.g.haps.f <- hcut(g.haps.f$dmat, k=k)
res.g.haps.m <- hcut(g.haps.m$dmat, k=k)

# dendrograms
lab_adj <- 1 # smaller -> more room for labels
cex <- 0.6
k_colors <- NULL
# uncorrected data
d1.1 <- fviz_dend(res.u.primates.f, k_colors=k_colors, labels_track_height=max(res.u.primates.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Primate females")
d1.2 <- fviz_dend(res.u.primates.m, k_colors=k_colors, labels_track_height=max(res.u.primates.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Primate males")
d1.3 <- fviz_dend(res.u.streps.f, k_colors=k_colors, labels_track_height=max(res.u.streps.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Strepsirrhine females")
d1.4 <- fviz_dend(res.u.streps.m, k_colors=k_colors, labels_track_height=max(res.u.streps.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Strepsirrhine males")
d1.5 <- fviz_dend(res.u.haps.f, k_colors=k_colors, labels_track_height=max(res.u.haps.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Haplorhine females")
d1.6 <- fviz_dend(res.u.haps.m, k_colors=k_colors, labels_track_height=max(res.u.haps.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Haplorhine males")
grid.arrange(d1.1, d1.2, d1.3, d1.4, d1.5, d1.6, ncol=3, nrow=2, as.table=FALSE)
# corrected data
d2.1 <- fviz_dend(res.c.primates.f, k_colors=k_colors, labels_track_height=max(res.c.primates.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Primate females")
d2.2 <- fviz_dend(res.c.primates.m, k_colors=k_colors, labels_track_height=max(res.c.primates.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Primate males")
d2.3 <- fviz_dend(res.c.streps.f, k_colors=k_colors, labels_track_height=max(res.c.streps.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Strepsirrhine females")
d2.4 <- fviz_dend(res.c.streps.m, k_colors=k_colors, labels_track_height=max(res.c.streps.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Strepsirrhine males")
d2.5 <- fviz_dend(res.c.haps.f, k_colors=k_colors, labels_track_height=max(res.c.haps.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Haplorhine females")
d2.6 <- fviz_dend(res.c.haps.m, k_colors=k_colors, labels_track_height=max(res.c.haps.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Haplorhine males")
grid.arrange(d2.1, d2.2, d2.3, d2.4, d2.5, d2.6, ncol=3, nrow=2, as.table=FALSE)
# Goswami's method
d3.1 <- fviz_dend(res.g.primates.f, k_colors=k_colors, labels_track_height=max(res.g.primates.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Primate females")
d3.2 <- fviz_dend(res.g.primates.m, k_colors=k_colors, labels_track_height=max(res.g.primates.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Primate males")
d3.3 <- fviz_dend(res.g.streps.f, k_colors=k_colors, labels_track_height=max(res.g.streps.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Strepsirrhine females")
d3.4 <- fviz_dend(res.g.streps.m, k_colors=k_colors, labels_track_height=max(res.g.streps.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Strepsirrhine males")
d3.5 <- fviz_dend(res.g.haps.f, k_colors=k_colors, labels_track_height=max(res.g.haps.f$height)/lab_adj, horiz=TRUE, cex=cex, main="Haplorhine females")
d3.6 <- fviz_dend(res.g.haps.m, k_colors=k_colors, labels_track_height=max(res.g.haps.m$height)/lab_adj, horiz=TRUE, cex=cex, main="Haplorhine males")
grid.arrange(d3.1, d3.2, d3.3, d3.4, d3.5, d3.6, ncol=3, nrow=2, as.table=FALSE)

# optimum number of modules based on the gap statistic
K.max <- 10
# no control variables
with(clusGap(u.primates.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(u.primates.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(u.streps.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(u.streps.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(u.haps.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(u.haps.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
# include control variables
with(clusGap(c.primates.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(c.primates.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(c.streps.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(c.streps.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(c.haps.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(c.haps.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
# Goswami covariance, PPA
with(clusGap(g.primates.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(g.primates.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(g.streps.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 3
with(clusGap(g.streps.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 7
with(clusGap(g.haps.f$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1
with(clusGap(g.haps.m$mds, hcut, K.max), maxSE(Tab[,"gap"],Tab[,"SE.sim"])) # 1

# Goswami modularity tests on Ward's cluster results for k=6
Goswami_modularity_test(g.primates.f$lcm, k=6) # p < 0.001 ***
Goswami_modularity_test(g.primates.m$lcm, k=6) # p < 0.001 ***
Goswami_modularity_test(g.streps.f$lcm, k=6) # p < 0.001 ***
Goswami_modularity_test(g.streps.m$lcm, k=6) # p < 0.001 ***
Goswami_modularity_test(g.haps.f$lcm, k=6) # p < 0.001 ***
Goswami_modularity_test(g.haps.m$lcm, k=6) # p < 0.001 ***

########
# END ##
########
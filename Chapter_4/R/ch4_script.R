################
# PREPARATIONS #
################

# load packages
library("geomorph")
library("plyr")
library("tidyverse")
library("abind") 
library("ggtree")

# load custom functions
source("Chapter_4/R/functions/read.pp.R")
source("Chapter_4/R/functions/plot_wireframes.R")
source("Chapter_4/R/functions/plotGMPhyloMorphoSpace_axisflip.R")

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

# add centroid size to data
data.f$Csize <- gpa.f$Csize
data.m$Csize <- gpa.m$Csize

# create geomorph dataframe for males and females
gdf.f <- geomorph.data.frame(gpa.f, phy=tree.f, 
                             PhysicalAlignment=data.f$PhysicalAlignment,
                             BodyMass=data.f$BodyMass,
                             DimorphismIndex=data.f$DimorphismIndex,
                             Nocturnal=data.f$Nocturnal,
                             Gouging=data.f$Gouging,
                             Gouging2=data.f$Gouging2,
                             Folivore=data.f$Folivore)
gdf.m <- geomorph.data.frame(gpa.m, phy=tree.m, 
                             PhysicalAlignment=data.m$PhysicalAlignment,
                             BodyMass=data.m$BodyMass,
                             DimorphismIndex=data.m$DimorphismIndex,
                             Nocturnal=data.m$Nocturnal,
                             Gouging=data.m$Gouging,
                             Gouging2=data.m$Gouging2,
                             Folivore=data.m$Folivore)

##############################
# EFFECT OF ALIGNMENT METHOD #
##############################

# shape vs alignment method
procD.pgls(coords ~ PhysicalAlignment, phy=phy, data=gdf.f) # p = 0.663
procD.pgls(coords ~ PhysicalAlignment, phy=phy, data=gdf.m) # p = 0.977

#######################
# PHYLOGENETIC SIGNAL #
#######################

# multivariate Blomberg's K
physignal(landmarks.f, tree.f) # K = 0.47, p = 0.003
physignal(landmarks.m, tree.m) # K = 0.511, p = 0.003

# after controlling for centroid size
single.csize.f <- procD.pgls(coords ~ log(Csize), phy=phy, data=gdf.f) # 0.018 *
single.csize.m <- procD.pgls(coords ~ log(Csize), phy=phy, data=gdf.m) # 0.002 **
physignal(single.csize.f$pgls.residuals, tree.f) # K = 0.43, p = 0.001
physignal(single.csize.m$pgls.residuals, tree.m) # K = 0.41, p = 0.001

# after controlling for all predictors
multi.f <- procD.pgls(coords ~  Gouging + Folivore + Nocturnal + DimorphismIndex + log(Csize), phy=phy, data=gdf.f) # p = 0.026 * **
multi.m <-procD.pgls(coords ~  Gouging + Folivore + Nocturnal + DimorphismIndex+ log(Csize), phy=phy, data=gdf.m) # p = 0.089 *
physignal(multi.f$pgls.residuals, tree.f) # K = 0.33, p = 0.001
physignal(multi.m$pgls.residuals, tree.m) # K = 0.34, p = 0.001

#####################
# INTERACTION TERMS #
#####################

# females
procD.pgls(coords ~ log(Csize) + DimorphismIndex*Nocturnal, phy=phy, data=gdf.f) # p = 0.467
procD.pgls(coords ~ log(Csize) + DimorphismIndex*Folivore, phy=phy, data=gdf.f) # p = 0.413
procD.pgls(coords ~ log(Csize) + DimorphismIndex*Gouging, phy=phy, data=gdf.f) # p = 0.653
procD.pgls(coords ~ log(Csize) + Gouging*Folivore, phy=phy, data=gdf.f) # fails
procD.pgls(coords ~ log(Csize) + Gouging2*Folivore, phy=phy, data=gdf.f) # fails
procD.pgls(coords ~ log(Csize) + Gouging*Nocturnal, phy=phy, data=gdf.f) # p = 0.581
procD.pgls(coords ~ log(Csize) + Folivore*Nocturnal, phy=phy, data=gdf.f) # p = 0.131

# males
procD.pgls(coords ~ log(Csize) + DimorphismIndex*Nocturnal, phy=phy, data=gdf.m) # p = 0.322
procD.pgls(coords ~ log(Csize) + DimorphismIndex*Folivore, phy=phy, data=gdf.m) # p = 0.814
procD.pgls(coords ~ log(Csize) + DimorphismIndex*Gouging, phy=phy, data=gdf.m) # p = 0.705
procD.pgls(coords ~ log(Csize) + Gouging*Folivore, phy=phy, data=gdf.m) # fails
procD.pgls(coords ~ log(Csize) + Gouging2*Folivore, phy=phy, data=gdf.m) # fails
procD.pgls(coords ~ log(Csize) + Gouging*Nocturnal, phy=phy, data=gdf.m) # p = 0.864
procD.pgls(coords ~ log(Csize) + Folivore*Nocturnal, phy=phy, data=gdf.m) # p = 0.437

##################################
# HYPOTHESIS TESTS - ALL SPECIES #
##################################

# single variables: females
single.noc.f <- procD.pgls(coords ~ log(Csize) + Nocturnal, phy=phy, data=gdf.f) # 0.292
single.fol.f <- procD.pgls(coords ~ log(Csize) + Folivore, phy=phy, data=gdf.f) # 0.036 *
single.gouge.f <- procD.pgls(coords ~ log(Csize) + Gouging, phy=phy, data=gdf.f) # 0.378
single.gouge.f2 <- procD.pgls(coords ~ log(Csize) + Gouging2, phy=phy, data=gdf.f) # 0.552
single.di.f <- procD.pgls(coords ~ log(Csize) + DimorphismIndex, phy=phy, data=gdf.f) # 0.001 **

# single variables: males
single.noc.m <- procD.pgls(coords ~ log(Csize) + Nocturnal, phy=phy, data=gdf.m) # 0.314
single.fol.m <- procD.pgls(coords ~ log(Csize) + Folivore, phy=phy, data=gdf.m) # 0.039 *
single.gouge.m <- procD.pgls(coords ~ log(Csize) + Gouging, phy=phy, data=gdf.m) # 0.746
single.gouge.m2 <- procD.pgls(coords ~ log(Csize) + Gouging2, phy=phy, data=gdf.m) # 0.661
single.di.m <- procD.pgls(coords ~ log(Csize) + DimorphismIndex, phy=phy, data=gdf.m) # 0.001 **

# full model: females
multi.csize.f <- procD.pgls(coords ~ DimorphismIndex + Gouging + Folivore + Nocturnal + log(Csize), phy=phy, data=gdf.f) # p = 0.026 *
multi.noc.f <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Gouging + Folivore + Nocturnal, phy=phy, data=gdf.f) # p = 0.241
multi.fol.f <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Gouging + Nocturnal + Folivore, phy=phy, data=gdf.f) # p = 0.068 .
multi.gouge.f <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Folivore + Nocturnal + Gouging, phy=phy, data=gdf.f) # p = 0.268
multi.gouge.f2 <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Folivore + Nocturnal + Gouging2, phy=phy, data=gdf.f) # p = 0.465
multi.di.f <- procD.pgls(coords ~ log(Csize) + Gouging + Folivore + Nocturnal + DimorphismIndex, phy=phy, data=gdf.f) # p = 0.001 **

# full model: males
multi.csize.m <- procD.pgls(coords ~ DimorphismIndex + Gouging + Folivore + Nocturnal + log(Csize), phy=phy, data=gdf.m) # p = 0.089 .
multi.noc.m <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Gouging + Folivore + Nocturnal, phy=phy, data=gdf.m) # p = 0.204
multi.fol.m <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Gouging + Nocturnal + Folivore, phy=phy, data=gdf.m) # p = 0.055 .
multi.gouge.m <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Folivore + Nocturnal + Gouging, phy=phy, data=gdf.m) # p = 0.676
multi.gouge.m2 <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Folivore + Nocturnal + Gouging2, phy=phy, data=gdf.m) # p = 0.687
multi.di.m <- procD.pgls(coords ~ log(Csize) + Gouging + Folivore + Nocturnal + DimorphismIndex, phy=phy, data=gdf.m) # p = 0.001 **

#########################################
# HYPOTHESIS TESTS - SUBSETS OF SPECIES #
#########################################

# effect of activity pattern below 75 mm skull length threshold (Kay & Cartmill 1977)
# compute skull lengths for females and identify species to be included
skull.lengths <- sqrt(colSums((landmarks.f["Glabella",,] - landmarks.f["Opisthocranium",,])^2))
subset.f.75 <- subset.m.75 <- names(skull.lengths[skull.lengths < 75]) # 44 taxa 
subset.m.75[subset.m.75 == "Cephalopachus_bancanus"] <- "Carlito_syrichta"
# subset and align shape, subset data and tree, create geomorph dataframe
gpa.f.75 <- gpagen(landmarks.f[,,subset.f.75])
gpa.m.75 <- gpagen(landmarks.m[,,subset.m.75])
data.f.75 <- filter(data.f, genus_species %in% subset.f.75)
data.m.75 <- filter(data.m, genus_species %in% subset.m.75)
tree.f.75 <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% subset.f.75)])
tree.m.75 <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% subset.m.75)])
# create geomorph object
gdf.f.75 <- geomorph.data.frame(gpa.f.75, phy=tree.f.75, 
                                PhysicalAlignment=data.f.75$PhysicalAlignment,
                                BodyMass=data.f.75$BodyMass,
                                DimorphismIndex=data.f.75$DimorphismIndex,
                                Nocturnal=data.f.75$Nocturnal,
                                Gouging=data.f.75$Gouging,
                                Folivore=data.f.75$Folivore)
gdf.m.75 <- geomorph.data.frame(gpa.m.75, phy=tree.m.75,
                                PhysicalAlignment=data.m.75$PhysicalAlignment,
                                BodyMass=data.m.75$BodyMass,
                                DimorphismIndex=data.m.75$DimorphismIndex,
                                Nocturnal=data.m.75$Nocturnal,
                                Gouging=data.m.75$Gouging,
                                Folivore=data.m.75$Folivore)
# simple model
single.noc.f75 <- procD.pgls(coords ~ log(Csize) + Nocturnal, phy=phy, data=gdf.f.75) # 0.212
single.noc.m75 <- procD.pgls(coords ~ log(Csize) + Nocturnal, phy=phy, data=gdf.m.75) # 0.263
# full model
multi.noc.f75 <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Gouging + Folivore + Nocturnal, phy=phy, data=gdf.f.75) # p = 0.173
multi.noc.m75 <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Gouging + Folivore + Nocturnal, phy=phy, data=gdf.m.75) # p = 0.139


# effect of folivory above Kay's threshold (500 g) (Kay 1984)
# subset to species >500 g
subset.f.500 <- subset.m.500 <- unlist(filter(data.f, BodyMass > 0.5) %>% select(genus_species)) # 50 taxa
subset.m.500 <- c(subset.f.500, "Rhinopithecus_roxellana") # 51 taxa
# subset and align shape, subset data and tree, create geomorph dataframe
gpa.f.500 <- gpagen(landmarks.f[,,subset.f.500])
gpa.m.500 <- gpagen(landmarks.m[,,subset.m.500])
data.f.500 <- filter(data.f, genus_species %in% subset.f.500)
data.m.500 <- filter(data.m, genus_species %in% subset.m.500)
tree.f.500 <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% subset.f.500)])
tree.m.500 <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% subset.m.500)])
# create geomorph object
gdf.f.500 <- geomorph.data.frame(gpa.f.500, phy=tree.f.500, 
                                PhysicalAlignment=data.f.500$PhysicalAlignment,
                                BodyMass=data.f.500$BodyMass,
                                DimorphismIndex=data.f.500$DimorphismIndex,
                                Nocturnal=data.f.500$Nocturnal,
                                Gouging=data.f.500$Gouging,
                                Folivore=data.f.500$Folivore)
gdf.m.500 <- geomorph.data.frame(gpa.m.500, phy=tree.m.500,
                                PhysicalAlignment=data.m.500$PhysicalAlignment,
                                BodyMass=data.m.500$BodyMass,
                                DimorphismIndex=data.m.500$DimorphismIndex,
                                Nocturnal=data.m.500$Nocturnal,
                                Gouging=data.m.500$Gouging,
                                Folivore=data.m.500$Folivore)
# simple model
single.fol.f500 <- procD.pgls(coords ~ log(Csize) + Folivore, phy=phy, data=gdf.f.500) # p = 0.043 *
single.fol.m500 <- procD.pgls(coords ~ log(Csize) + Folivore, phy=phy, data=gdf.m.500) # p = 0.043 *
# full model
multi.fol.f500 <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Gouging + Nocturnal + Folivore, phy=phy, data=gdf.f.500) # p = 0.107
multi.fol.m500 <- procD.pgls(coords ~ log(Csize) + DimorphismIndex + Gouging + Nocturnal + Folivore, phy=phy, data=gdf.m.500) # p = 0.039 *


############
# 2D PLOTS #
############

# log centroid size vs log body mass
CS_BM <- data.frame(LogCentroid = log(c(gpa.f$Csize, gpa.m$Csize)), 
                    LogBody = log(c(data.f$BodyMass, data.m$BodyMass)),
                    Sex = c(rep("Females",nrow(data.f)),rep("Males",nrow(data.m))),
                    Clade = c(as.character(data.f$Superfamily), as.character(data.m$Superfamily)))
ggplot(CS_BM, aes(x=LogBody, y=LogCentroid, fill=Clade, color=Clade)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha=0.1) +
  facet_grid(~Sex) +
  ylab("Log centroid size") +
  xlab("Log body mass") + 
  theme(axis.text.x=element_text(size = 12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(20,0,0,0), size = 14, face = 2),
        axis.title.y=element_text(margin=margin(0,20,0,0), size = 14, face = 2),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        strip.text=element_text(size=14))

# sexual dimorphism vs log body mass
DI_BM <- ddply(data, .(genus_species), function(x) {
  data.frame(LogDi = x$DimorphismIndex[1], 
             LogBody = log(x$BodyMass)[1],
             Clade = x$Superfamily[1])
})
ggplot(DI_BM, aes(x=LogBody, y=LogDi, fill=Clade, color=Clade)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha=0.1) +
  ylab("Dimorphism index") +
  xlab("Log body mass") + 
  theme(axis.text.x=element_text(size = 12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(20,0,0,0), size = 14, face = 2),
        axis.title.y=element_text(margin=margin(0,20,0,0), size = 14, face = 2),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

# phylogeny with ecological / behavioral data
treefig.dat <- ddply(data, .(genus_species), function (x) {
  temp <- x[1,c("Nocturnal","Folivore","Gouging","Gouging2")]
  ifelse(temp == 0, "Absent", "Present")
})
rownames(treefig.dat) <- sub("_"," ", treefig.dat$genus_species)
treefig.dat <- treefig.dat[,2:5]
tree2 <- tree
tree2$tip.label <- sub("_"," ",tree2$tip.label)
fig <- ggplot(tree2) + 
  geom_tree() + 
  theme_tree() + 
  geom_tiplab(list(tree2$tip.label), size=2.5, offset=7, fontface="italic") + 
  xlim(0,95) +
  ylim(0,80)
fig <- gheatmap(fig, treefig.dat, width=0.08, colnames_position="top", colnames_angle=85, colnames_offset_y=0.5, hjust=0, font.size=3, color="black") +
  scale_fill_manual(breaks=c("Absent","Present"), values=c("white", "black"))
fig <- phylopic(fig, "f598fb39-facf-43ea-a576-1861304b2fe4", node=110) # tarsiers
fig <- phylopic(fig, "aceb287d-84cf-46f1-868c-4797c4ac54a8", node=96) # nwm
fig <- phylopic(fig, "bac25f49-97a4-4aec-beb6-f542158ebd23", node=113) # lemurs
fig <- phylopic(fig, "7fb9bea8-e758-4986-afb2-95a2c3bf983d", node=124) # galagoids
fig <- phylopic(fig, "72f2f854-f3cd-4666-887c-35d5c256ab0f", node=70) # owm
fig <- phylopic(fig, "0174801d-15a6-4668-bfe0-4c421fbe51e8", node=89) # hominoids
fig + theme(legend.position=c(0.15,0.85), legend.text=element_text(size=13))

# phylomorphospace
layout(matrix(1:2, 2, 1))
plotGMPhyloMorphoSpace_axisflip(tree.f, gpa.f$coords, tip.text=gsub("_.*","",tree.f$tip.label), node.labels=F, plot.param=list(t.cex=0.3, n.cex=0.3, lwd=0.3, txt.cex=0.8))
mtext("Females", line=1, cex=1.5)
plotGMPhyloMorphoSpace_axisflip(tree.m, gpa.m$coords, tip.text=gsub("_.*","",tree.m$tip.label), node.labels=F, yaxis=-2, plot.param=list(t.cex=0.3, n.cex=0.3, lwd=0.3, txt.cex=0.8))
mtext("Males", line=1, cex=1.5)

######################
# 3D WIREFRAME PLOTS #
######################

# load wireframe links
wireframe <- read.csv("./Chapter_4/data/wireframe.csv", header=FALSE)
lines.col <- mapvalues(wireframe[,1], unique(wireframe[,1]), c("red","green","blue","purple","goldenrod"))

# mean shape
layout3d(matrix(1:2, 1, 2), sharedMouse = TRUE)
plot.coords(gpa.f$consensus, wireframe[,2:3], lines.col=lines.col, add=TRUE, 
            legend=c("Mandible","Face","Braincase","Zygomatic","Basicranium"), legend.pos="topright", 
            legend.col=c("red","green","blue","purple","goldenrod"))
plot.coords(gpa.m$consensus, wireframe[,2:3], lines.col=lines.col)
snapshot3d(filename = '~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Correlates of skull shape/Figures/Figure 5.5 - mean skull shape.png', fmt = 'png')

# centroid size
layout3d(matrix(1:4, 2, 2, byrow=TRUE), sharedMouse = TRUE)
plot.procD(multi.csize.f, wireframe[,2:3], value=log(min(gdf.f$Csize)))
plot.procD(multi.csize.f, wireframe[,2:3], value=log(max(gdf.f$Csize)))
plot.procD(multi.csize.m, wireframe[,2:3], value=log(min(gdf.m$Csize)))
plot.procD(multi.csize.m, wireframe[,2:3], value=log(max(gdf.m$Csize)))
snapshot3d(filename = '~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Correlates of skull shape/Figures/Figure 5.6 - centroid size wireframes.png', fmt = 'png')

# dimorphism 
layout3d(matrix(1:4, 2, 2, byrow=TRUE), sharedMouse = TRUE)
plot.procD(multi.di.f, wireframe[,2:3], value=min(gdf.f$DimorphismIndex))
plot.procD(multi.di.f, wireframe[,2:3], value=max(gdf.f$DimorphismIndex))
plot.procD(multi.di.m, wireframe[,2:3], value=min(gdf.m$DimorphismIndex))
plot.procD(multi.di.m, wireframe[,2:3], value=max(gdf.m$DimorphismIndex))
snapshot3d(filename = '~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Correlates of skull shape/Figures/Figure 5.7 - dimorphism wireframes.png', fmt = 'png')

# nocturnality 
layout3d(matrix(1:2, 1, 2), sharedMouse = TRUE)
plot.procD(multi.noc.f, wireframe[,2:3], value=0, points.col="black", lines.col="black", 
           legend=c("Non-nocturnal","Nocturnal"), legend.pos="topright", legend.col=c("black","skyblue"))
plot.procD(multi.noc.f, wireframe[,2:3], value=1, points.col="skyblue", lines.col="skyblue", add=TRUE)
plot.procD(multi.noc.m, wireframe[,2:3], value=0, points.col="black", lines.col="black")
plot.procD(multi.noc.m, wireframe[,2:3], value=1, points.col="skyblue", lines.col="skyblue", add=TRUE)
snapshot3d(filename = '~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Correlates of skull shape/Figures/Figure 5.8 - nocturnality wireframes.png', fmt = 'png')

# folivory
layout3d(matrix(1:2, 1, 2), sharedMouse = TRUE)
plot.procD(multi.fol.f, wireframe[,2:3], value=0, points.col="black", lines.col="black", 
           legend=c("Non-folivorous","Folivorous"), legend.pos="topright", legend.col=c("black","palegreen"))
plot.procD(multi.fol.f, wireframe[,2:3], value=1, points.col="palegreen", lines.col="palegreen", add=TRUE)
plot.procD(multi.fol.m, wireframe[,2:3], value=0, points.col="black", lines.col="black")
plot.procD(multi.fol.m, wireframe[,2:3], value=1, points.col="palegreen", lines.col="palegreen", add=TRUE)
snapshot3d(filename = '~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Correlates of skull shape/Figures/Figure 5.9 - folivory wireframes.png', fmt = 'png')

# gouging 
layout3d(matrix(1:2, 1, 2), sharedMouse = TRUE)
plot.procD(multi.gouge.f, wireframe[,2:3], value=0, points.col="black", lines.col="black", 
           legend=c("Non-gouging","Gouging"), legend.pos="topright", legend.col=c("black","pink"))
plot.procD(multi.gouge.f, wireframe[,2:3], value=1, points.col="pink", lines.col="pink", add=TRUE)
plot.procD(multi.gouge.m, wireframe[,2:3], value=0, points.col="black", lines.col="black")
plot.procD(multi.gouge.m, wireframe[,2:3], value=1, points.col="pink", lines.col="pink", add=TRUE)
snapshot3d(filename = '~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Correlates of skull shape/Figures/Figure 5.10 - gouging wireframes.png', fmt = 'png')

# gouging2
layout3d(matrix(1:2, 1, 2), sharedMouse = TRUE)
plot.procD(multi.gouge.f2, wireframe[,2:3], value=0, points.col="black", lines.col="black", 
           legend=c("Non-gouging","Gouging"), legend.pos="topright", legend.col=c("black","pink"))
plot.procD(multi.gouge.f2, wireframe[,2:3], value=1, points.col="pink", lines.col="pink", add=TRUE)
plot.procD(multi.gouge.m2, wireframe[,2:3], value=0, points.col="black", lines.col="black")
plot.procD(multi.gouge.m2, wireframe[,2:3], value=1, points.col="pink", lines.col="pink", add=TRUE)
snapshot3d(filename = '~/Dropbox/Dissertation/Dissertation chapters/Chapter 4 - Correlates of skull shape/Figures/Figure 5.10 - gouging2 wireframes.png', fmt = 'png')

# nocturnality in species <75 mm skull length
layout3d(matrix(1:2, 1, 2), sharedMouse = TRUE)
plot.procD(single.noc.f75, wireframe[,2:3], value=0, points.col="black", lines.col="black", 
           legend=c("Non-nocturnal","Nocturnal"), legend.pos="topright", legend.col=c("black","skyblue"))
plot.procD(single.noc.f75, wireframe[,2:3], value=1, points.col="skyblue", lines.col="skyblue", add=TRUE)
plot.procD(single.noc.m75, wireframe[,2:3], value=0, points.col="black", lines.col="black")
plot.procD(single.noc.m75, wireframe[,2:3], value=1, points.col="skyblue", lines.col="skyblue", add=TRUE)

# folivory in species >500g body mass
layout3d(matrix(1:2, 1, 2), sharedMouse = TRUE)
plot.procD(single.fol.f500, wireframe[,2:3], value=0, points.col="black", lines.col="black", 
           legend=c("Non-folivorous","Folivorous"), legend.pos="topright", legend.col=c("black","palegreen"))
plot.procD(single.fol.f500, wireframe[,2:3], value=1, points.col="palegreen", lines.col="palegreen", add=TRUE)
plot.procD(single.fol.m500, wireframe[,2:3], value=0, points.col="black", lines.col="black")
plot.procD(single.fol.m500, wireframe[,2:3], value=1, points.col="palegreen", lines.col="palegreen", add=TRUE)

########
# END ##
########
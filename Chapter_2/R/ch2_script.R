################
# PREPARATIONS #
################

# load packages
library("plyr")
library("ape")
library("ggplot2")
library("ggtree")
library("EBImage")

# read data and tree
data <- read.csv("Chapter_2/data/table1.csv")
tree <- read.nexus("Chapter_2/data/tree.nex")

#################
# SUMMARY STATS #
#################

# number of rows 
nrow(data)

# number of known males/females
table(data$Sex)

# proportion of target sample size (N = 72) attained
round((sum(data$SexDat %in% c("M", "B"))/72)*100, 1) # proportion of males
round((sum(data$SexDat %in% c("F", "B"))/72)*100, 1) # proportion of females

# scan sources
table(data$Scan.source)
round((sum(table(data$Scan.source)[2:4])/nrow(data)), 2) # proportion publicly availble
round((table(data$Scan.source)[1]/nrow(data)), 2) # proportion de novo
round((table(data$Scan.source)[5]/nrow(data)), 2) # proportion personal communication

# sources of publicly available scans
table(data$Scan.source)[2]/sum(table(data$Scan.source)[2:4]) # KUPRI
table(data$Scan.source)[3]/sum(table(data$Scan.source)[2:4]) # smithsonian 3D collection
table(data$Scan.source)[4]/sum(table(data$Scan.source)[2:4]) # morphosource

####################
# PHYLOGENY FIGURE #
####################

# check that names all match up
test <- paste(data$Genus, data$Species, sep="_")
test[!test %in% tree$tip.label]
tree$tip.label[!tree$tip.label %in% test]
rm(test)

# set tree labels to genus name only
tree$tip.label <- sapply(strsplit(tree$tip.label, "_"), "[[", 1)

# create table for heatmap
tab <- ddply(data, .(Genus), function(x) {
  m <- ifelse("M" %in% x$SexDat | "B" %in% x$SexDat, 1, 0)
  f <- ifelse("F" %in% x$SexDat | "B" %in% x$SexDat, 1, 0)
  data.frame(M=m, F=f)
})
row.names(tab) <- tab$Genus
tab <- tab[,-1]

# plot tree
fig <- ggplot(tree) + 
  geom_tree() + 
  theme_tree() + 
  geom_tiplab(list(tree$tip.label), size=2.5, offset=5) + 
  xlim(0,95) +
  geom_cladelabel(96, "Ceboidea", offset= 18, barsize=2, angle=90, offset.text = 1.2, hjust=0.5) +
  geom_cladelabel(70, "Cercopithecoidea", offset= 18, barsize=2, angle=90, offset.text = 1, hjust=0.5) +
  geom_cladelabel(113, "Lemuroidea", offset= 18, barsize=2, angle=90, offset.text = 1.2, hjust=0.5) +
  geom_cladelabel(124, "Galagoidea", offset= 18, barsize=2, angle=90, offset.text = 1.2, hjust=0.5) +
  geom_cladelabel(89, "Hominoidea", offset= 18, barsize=2, angle=90, offset.text = 1.2, hjust=0.5) +
  geom_cladelabel(110, "Tarsioidea", offset= 18, barsize=2, angle=90, offset.text = 1.2, hjust=0.5, fontsize = 3)
fig <- gheatmap(fig, tab, offset=0.1, width=0.05, low="white", high="black", colnames_position = "top", font.size=2.5)
fig <- phylopic(fig, "f598fb39-facf-43ea-a576-1861304b2fe4", node=110) # tarsiers
fig <- phylopic(fig, "aceb287d-84cf-46f1-868c-4797c4ac54a8", node=96) # nwm
fig <- phylopic(fig, "bac25f49-97a4-4aec-beb6-f542158ebd23", node=113) # lemurs
fig <- phylopic(fig, "7fb9bea8-e758-4986-afb2-95a2c3bf983d", node=124) # galagoids
fig <- phylopic(fig, "72f2f854-f3cd-4666-887c-35d5c256ab0f", node=70) # owm
fig <- phylopic(fig, "0174801d-15a6-4668-bfe0-4c421fbe51e8", node=89) # hominoids
fig + theme(legend.position="none")

########
# END ##
########
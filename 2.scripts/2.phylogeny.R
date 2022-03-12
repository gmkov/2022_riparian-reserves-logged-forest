################################################################################################
### Turning a genus-level phylogeny into a species-level phylogeny with genus-level polytomies ###
###   fig 4 tree prep ###
################################################################################################
# load packages
rm(list=ls())
dev.off()
library("ape")
library("Biostrings")
library("ggtree")
library("ggplot2")
library("ggtree")
library(tidytree)
library(treeio)
library("phytools")
library("adephylo")
library("phylobase")
library(rotl)
library(httr)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)

#BiocManager::install("ggtree")
#install.packages("remotes")
#install.packages("XML", repos = "http://www.omegahat.net/R")
#########################  0. functions ########### 
# from https://github.com/willpearse/phyloGenerator/blob/master/merging.R
require(ape)
#Replace a given list of genera in a phylogeny with polytomes of their species
#Needs a tree, and vectors of genera and species
#NOTE: vectors should be of same length, e.g. c('Quercus', 'Quercus', 'Phocoena') and c('ilex', 'robur', 'phocoena')


make.composite.with.polytomies <- function(tree, genera, species, max.genus.age=NA){
  #Functions#
  #Binds a clade into a phylogeny, replacing a given tip, and making it ultrametric	
  #A hack to fix the ultrametric nature, but APE is a pile of shite so it's all I can do	
  bind.ultrametric.by.replacement <- function(backbone, donor, replacing.tip.label, donor.length=NA){	
    #Find the species we're replacing	
    bind.point <- which(backbone$tip.label == replacing.tip.label)	
    #Bind *badly* the two trees together	
    backbone <- bind.tree(backbone, donor, where=bind.point)	
    #Now we've bound, where's the genus located?	
    which.tip <- which(backbone$tip.label == donor$tip.label[1])	
    #What node is above the genus?	
    which.node <- backbone$edge[which(backbone$edge[,2] == which.tip),1]	
    #Where is that node in the edge.list (and thus edge.length list)	
    which.edge <- which(backbone$edge[,2] == which.node)	
    #What length is that branch?	
    tip.length <- backbone$edge.length[which.edge]
    if(is.na(donor.length)){
      #It's twice as long as it should be, so replace it
      backbone$edge.length[which.edge] <- tip.length/2
    } else {
      #Correct for the manual adjustment...
      backbone$edge.length[which.edge] <- tip.length - donor.length/2
    }
    #Return the ultrametric tree!	
    return(backbone)	
  }	
  #Make a polytomy from a given species list, optionally with a tip.length	
  make.polytomy <- function(species, tip.length=NA){	
    d.f <- data.frame(spp=factor(species))	
    polytomy <- as.phylo.formula(~spp, data=d.f)	
    if(!is.na(tip.length)) polytomy$edge.length <- rep(tip.length, length(species))	
    return(polytomy)	
  }	
  #Find the unique branch length of a tip on a given phlogeny	
  find.unique.branch.length <- function(tree, tip){	
    #What number of tip.label is the tip	
    which.tip <- which(tree$tip.label == tip)	
    #What edge is the tip.label on	
    which.edge <- which(tree$edge[,2] == which.tip)	
    #What's the edge.length for that edge?	
    tip.length <- tree$edge.length[which.edge]	
    #Return that edge.length	
    return(tip.length)	
  }	
  
  #Code#
  #The genera and species vectors *must* be characters, so convert
  genera<-as.character(genera);species<-as.character(species)
  #Go along all the genera to be replaced with polytomes
  for(genus in unique(genera)){
    #Make a vector of the species that will go into this particular polytomy
    species.to.bind <- species[genera == genus]
    #If the length of that vector is 1...
    if(length(species.to.bind) == 1){
      #...we don't need to bind anything in, so just change the tip's name
      tree$tip.label[tree$tip.label == genus] <- species.to.bind
    } else {
      #Other, find the branch length leading to the particular genus
      tip.length <- find.unique.branch.length(tree, genus)
      #Don't warn about edge issues unless we need to think about it...
      edge.warning <- NA
      #Set the genus to the correct age
      if(!is.na(max.genus.age)){
        if(max.genus.age*2 < tip.length){
          tip.length <- min(tip.length, max.genus.age*2)
          edge.warning <- tip.length
        }
      }
      #Make the polytomy of species, with branch lengths appropriate for that genus
      polytomy <- make.polytomy(species.to.bind, (tip.length/2))
      #Bind that polytomy in, ultrametrially, to that phylogeny
      tree <- bind.ultrametric.by.replacement(tree, polytomy, genus, edge.warning)
    }
  }
  #Return the tree
  return(tree)
}

######################### 1. nymph tree  -data preprare ############
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/git/2022_riparian/")
traps <- read.csv("1.data/Points.csv", h = TRUE)
sp <- read.csv("1.data/Butterflies.csv", h = TRUE)

# 3 genes 
tree <- read.nexus("1.data/tree.3genes.representatives.tre")

ggtree(tree) + 
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, colour="grey", size=0)  +
  geom_tiplab(align=TRUE,  linetype='dashed', linesize=0,colour="black", size=3)  + xlim(-0, 2)

# all species found
unique(sp$genus); unique(sp$species)
sp$genus_sp <-  paste(sp$genus, sp$sp, sep = "_")
sp <- subset(sp, family=="Nymphalidae")
sp.for.tree <- data.frame(species=unique(sp$genus_sp))
sp.for.tree$genus <- unlist(lapply(str_split(sp.for.tree$species, "_"), head, 1)); head(sp.for.tree)
sp.for.tree <- subset(sp.for.tree, genus!=""&genus!="NA")
genera <- unique(sp.for.tree$genus)

# create list with genus in tree
sp.in.tree <- c(tree$tip.label)

genus.in.tree <- unlist(lapply(str_split(sp.in.tree, "_"), head, 1))

# genera that are missing at least one species
sp.not.in.tree <- subset(sp.for.tree, !(species %in% sp.in.tree)); unique(sp.not.in.tree$genus)
genus.missing.sp.in.tree <-  data.frame(genus=unique(sp.not.in.tree$genus))

sp.per.genus <- summarise(group_by(sp.for.tree, genus), n=n())
genus.missing.sp.in.tree$no.sp.total <- sp.per.genus$n[match(genus.missing.sp.in.tree$genus, sp.per.genus$genus)]
sum((genus.missing.sp.in.tree$no.sp.total))

####### from sp tree with duplicates, to sp tree without duplicates, then genus, then polytomies ####### 
# keep only one species per genus (in the tree), pick first that comes up (lowest random value)
tree.sp <- data.frame(species= sp.in.tree, genus=genus.in.tree  )

# keep only memphis glauce
tree.sp <-subset(tree.sp, (genus!="Memphis"|species=="Memphis_glauce"))

#subset to unique
tree.sp$random.value <- 1:nrow(tree.sp)
tree.sp.one <- tree.sp %>% 
  group_by(genus) %>% 
  filter(random.value== min(random.value)) %>% 
  distinct
tree.sp.one$species <- droplevels(tree.sp.one$species)

sp.to.drop <- subset(sp.in.tree, !(sp.in.tree%in%tree.sp.one$species ))

# subset species tree to match these, only one sp per genus
plot(drop.tip(tree, sp.to.drop))
tree.sp.one <- drop.tip(tree, sp.to.drop)

# make sp tree into genus tree
tree.genus <- tree.sp.one 
tree.genus$tip.label <- unlist(lapply(str_split(tree.genus$tip.label, "_"), head, 1))
plot(tree.genus)

# make genus tree into sp tree with polytomies, 116 species in tree/ 121
genus.sp.info1 <- subset(sp.for.tree, sp.for.tree$genus %in% tree.genus$tip.label , drop = T);genus.sp.info1

tree.genus.sp <- make.composite.with.polytomies(tree.genus, c(as.character(genus.sp.info1$genus)), 
                                    c(as.character(genus.sp.info1$species)) )

# drop tips that correspond to species in Peru, not in jamari. they didnt change to species name, so only genus
sp.to.drop.not.jamari <- subset(sp.in.tree, !(sp.in.tree %in% genus.sp.info1$species))
tree.genus.sp <- drop.tip(tree.genus.sp ,unlist(lapply(str_split(sp.to.drop.not.jamari, "_"), head, 1)) ); plot(tree.genus.sp)

unique(unlist(lapply(str_split(tree.genus.sp$tip.label, "_"), head, 1))) #50
unique(unlist(lapply(str_split(tree.genus.sp$tip.label, "_"), tail, 1))) #115
unique(sp$genus); unique(sp$species) # 52 genera, 144 species

unique(sp$genus)[!(unique(sp$genus) %in% unique(unlist(lapply(str_split(tree.genus.sp$tip.label, "_"), head, 1))))]
unique(sp$sp)[!(unique(sp$sp) %in% unique(unlist(lapply(str_split(tree.genus.sp$tip.label, "_"), tail, 1))))]

# write.nexus(tree.genus.sp, file="1.data/tree.3genes.genus.sp.nex")

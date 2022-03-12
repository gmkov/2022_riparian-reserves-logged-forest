########### RIPARIAN STRIPS analyses #######
rm(list=ls())
dev.off()

### PACKAGES ####
library(dplyr)
library(tidyr)
library(rgdal)
library(iNEXT)
library(ggplot2)
library(vegan)
library(ggpubr)
library(sf)
library(sp)
library(ggmap)
require(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)
library(CommEcol)
library(betapart)
library(cowplot)
library(phylosignal)
library(treeio)
library(phylobase)
library(adephylo)
library(ggbeeswarm)
library(stringr)
library(dplyr)

########################## DATA ##########
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/git/2022_riparian/")
traps <- read.csv("1.data/Points.csv", h = TRUE)
sp <- read.csv("1.data/Butterflies.csv", h = TRUE)

rivers <- readOGR("1.data/shapefiles",layer="rivers")
roads <- readOGR("1.data/shapefiles",layer="roads")
upa <- readOGR("1.data/shapefiles",layer="upa")
points <- readOGR("1.data/shapefiles",layer="points")
veg <- read.csv("1.data/Veg.csv", h = TRUE)

#### subsetting ####
# rm recaptures, and escaped
sp <- sp[sp$recaptured == "No", ]
sp <- sp[!is.na(sp$genus_sp), ]
sp <- sp[sp$genus_sp!="?_?", ]; sp <- sp[sp$genus_sp!="NA_NA", ]
length(unique(sp$genus_sp))
sp <- sp[sp$family=="Nymphalidae", ]

# matrices at two levels: strip (n=16) or point (n=64)
sp$habitat <-  paste(substr(sp$upa_id, 0, 3), substr(sp$location, 0,3), sep=".")
sp$ind <- 1; sp.strip <- as.data.frame(sp %>% group_by(upa_id, site_id, strip_id, stream, location,habitat , genus_sp) %>%
                                         dplyr::summarise(ind = sum(ind)) %>%
                                         spread(genus_sp, ind, fill = 0))

sp$ind <- 1;sp.point <- as.data.frame(sp %>% group_by(upa_id, site_id, strip_id, point_id, stream, location,habitat , genus_sp) %>%
                                        dplyr::summarise(ind = sum(ind)) %>%
                                        spread(genus_sp, ind, fill = 0)); sp.point

sp$ind <- 1;sp.habitat <- as.data.frame(sp %>% group_by(location,habitat , genus_sp) %>%
                                          dplyr::summarise(ind = sum(ind)) %>%
                                          spread(genus_sp, ind, fill = 0)); sp.habitat

# create a pair id
sp.point$pair.id <- paste(substr(sp.point$point_id,0,6),substr(sp.point$point_id,8,9), sep = "" ); sp.point$pair.id


#### abundance and species richness differences, summaries
require(dplyr)
# summ.hab.loc <-NA dplyr::summarise(dplyr::group_by(sp, upa_id, location,habitat),
#                                  abundance=n(),
#                                  sp.richness=length(unique(genus_sp)),
#                                  upa.id=unique(upa_id)); summ.hab.loc 
# 
# summ.hab.loc$habitat.2 <- ifelse(substr(summ.hab.loc$upa_id,0,1)=="U", "logged", "primary")

summ.point <- dplyr::summarise(group_by(sp, point_id, location, stream, upa_id, habitat),
                               abundance=n(),
                               sp.richness=length(unique(genus_sp))); summ.point
summ.point$habitat.2 <- ifelse(substr(summ.point$upa_id,0,1)=="U", "logged", "primary")

summ.strip<- dplyr::summarise(group_by(sp, strip_id, location, stream, upa_id, habitat),
                               abundance=n(),
                               sp.richness=length(unique(genus_sp))); summ.strip
summ.strip$habitat.2 <- ifelse(substr(summ.strip$upa_id,0,1)=="U", "logged", "primary")
summ.strip$habitat.3 <- ifelse(substr(summ.strip$habitat.2,0,1)=="p", "primary", summ.strip$habitat)
summ.point$habitat.3 <- ifelse(substr(summ.point$habitat.2,0,1)=="p", "primary", summ.point$habitat)

### veg data per point ####
veg.point <- dplyr::summarise(group_by(veg, point_code, site, location, stream),
                               max.canopy=mean(max_canopy),
                              mean.canopy=mean(mean_canopy),
                              canopy.cover=mean(canopy_cover),
                              tree.count=sum(tree_count),
                              diam.tree.total=sum(diam_sum),
                              diam.tree.mean=diam.tree.total/tree.count,
                              stream.width=mean(stream_width),
                              channel.width=mean(channel_width),
                              stream.depth=mean(stream_depth)); veg.point
summ.point$point_id==veg.point$point_code
veg.point$habitat.location <- paste(substr(veg.point$site,0,3), veg.point$location, sep=".")

### veg data per strip ####
veg$strip_id <- sp.point$strip_id[match(veg$point_code, sp.point$point_id)]
veg.strip <- dplyr::summarise(group_by(veg, strip_id, site, location, stream),
                              max.canopy=mean(max_canopy),
                              mean.canopy=mean(mean_canopy),
                              canopy.cover=mean(canopy_cover),
                              tree.count=sum(tree_count),
                              diam.tree.total=sum(diam_sum),
                              diam.tree.mean=diam.tree.total/tree.count,
                              stream.width=mean(stream_width),
                              channel.width=mean(channel_width),
                              stream.depth=mean(stream_depth)); veg.strip
summ.strip$strip_id==veg.strip$strip_id
head(summ.strip)
veg.strip$habitat.location <- paste(substr(veg.strip$site,0,3), veg.strip$location, sep=".")


##########################  MAPS ##############
# jpeg("plots/map.upa.river.points.jpeg", width = 1800, height = 1600, res = 250)
plot(upa, axes = TRUE, col = as.numeric(upa@data$UPA) + 1,
     xlim = c(494000, 510000), ylim = c(8960000, 8966000), asp = 1)
plot(upa, add = TRUE, col = rgb(1, 1, 1, 0.9))
plot(roads, add = TRUE, col = rgb(0.9,0,0, 0.6))
plot(rivers, add = TRUE, col = "blue")
plot(points, add = TRUE, cex = .75)

dev.off()

### SA insert
# plot extent -63.051,-9.419 : -62.909,-9.351
world <- ne_countries(scale = "medium", returnclass = "sf")
theme_set(theme_bw())
ggplot(data = world) +
  geom_sf(color="darkgrey") + ylab("Latitude") + xlab("Longitude")+
  coord_sf( xlim=c(-90,-40), ylim=c(-20,20), expand = FALSE)+
  theme(panel.grid.major = element_line(color = "white", linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "white"))+
  #scale_color_viridis(name="Altitude\n(m.a.s.l)")+
  scale_y_continuous(breaks=c(-10,0,10))+
  scale_x_continuous(breaks=c(-85, -70,  -60))+
  geom_rect(mapping = aes(xmin=-63.051, xmax=-62.909, ymin=-9.419, ymax=-9.351),fill=NA,color="black", alpha=.5)+
  geom_rect(mapping = aes(xmin=-65.051, xmax=-60.909, ymin=-10.419, ymax=-8.351),fill=NA,color="black", alpha=.5, lty=3)

# ggsave("plots/0.maps/SA.insert.png", width = 6, height = 4, dpi = 300)

########################## FIG. 2 - ABUNDANCE SP. RICHNESS POINT LEVEL ################
######## Fig. 2A. abundance all rip/int point-level ######
viridisLite::viridis(20)
# "#009E73","#0072B2", "#56B4E9", "#CC79A7"
# only compare same with same, interiors vs inter
# or only compare those that are n.s. in GLM n.s. PRI.Str - PRI.Int & UPA.Str - UPA.Int
summ.point$hab.loc <- paste(summ.point$location, summ.point$habitat.2)
summ.point$hab.loc <- factor(summ.point$hab.loc, levels = c("Interior primary", "Stream primary","Stream logged", "Interior logged"))
my_comparisons <- list(c("Interior primary", "Stream primary"),
                       c("Stream primary","Stream logged"),
                       c("Stream logged", "Interior logged"), 
                       c("Interior primary", "Interior logged"))
# get means for text
summarise(group_by(summ.point, habitat),
          ab=mean(abundance),
          sp=mean(sp.richness))

# FIG 2 A
fig.2a <- ggplot(summ.point, aes(x=hab.loc, y=abundance))+
  #geom_rect(aes(xmin = 1.5, xmax = 3.5, ymin = 15, ymax = 80),fill = "lightblue", alpha = 0.007)+
  stat_boxplot(geom ='errorbar', width=0.15) + 
  geom_boxplot(aes(fill=hab.loc),outlier.size = 1)+
  #scale_fill_manual(values=c("#009E73","#009E73",  "#56B4E9", "#56B4E9"))+
  scale_fill_manual(values=c("#009E73","#0072B2", "#56B4E9", "#CC79A7"))+
  #geom_text(size=5,color="black",aes(x=0.65,y=60, label="B1")) +
  theme_bw()+ylim(15,107)+
  theme(legend.position="none")+
  stat_compare_means(method="t.test", comparisons = my_comparisons, label = "p.signif", tip.length = 0.0,vjust = .1 , step.increase = 0.05)+
  labs(x="") +
  labs(y="Abundance")+
  #scale_x_discrete(labels=c("Interior\nprimary", "Riparian\nprimary" ,"Riparian\nreserve\nlogged", "Interior\nlogged"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        strip.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        axis.text =element_text(size = 12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.text.x = element_text(size=14, face = c("plain", "bold", "bold","plain")),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(colour = c("black","blue","blue","black"), size = c(1,7,7,1))); fig.2a

######## Fig. 2B. sp. richness all rip/int prim/log point-level ######
fig.2b <- ggplot(summ.point, aes(x=hab.loc, y=sp.richness))+
  #geom_rect(aes(xmin = 1.5, xmax = 3.5, ymin = 15, ymax = 80),fill = "lightblue", alpha = 0.007)+
  stat_boxplot(geom ='errorbar', width=0.15) + 
  geom_boxplot(aes(fill=hab.loc),outlier.size = 1)+
  scale_fill_manual(values=c("#009E73","#009E73",  "#56B4E9", "#56B4E9"))+
  scale_fill_manual(values=c("#009E73","#0072B2", "#56B4E9", "#CC79A7"))+
  #geom_text(size=5,color="black",aes(x=0.65,y=60, label="B1")) +
  theme_bw()+ylim(5,44)+
  theme(legend.position="none")+
  stat_compare_means(method="t.test", comparisons = my_comparisons,label = "p.signif", tip.length = 0.0,vjust = .1 , step.increase = 0.05)+
  labs(x="Forest") + 
  #ylim(8,45)+
  labs(y="Species richness")+
  scale_x_discrete(labels=c("Interior\nprimary", "Riparian\nprimary" ,"Riparian\nreserve\nlogged", "Interior\nlogged"))+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text.x = element_text(size = 16),
        axis.text =element_text(size = 12), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_text(size=14, face = c("plain", "bold", "bold","plain")),
        axis.ticks.x = element_line(colour = c("black","blue","blue","black"), size = c(1,7,7,1))); fig.2b

cowplot::plot_grid(fig.2a, fig.2b, ncol = 1,align = "v", rel_heights = c(1,1.4), labels = c("A","B"), label_size = 20)
#ggsave("plots/fig2.png", width = 4.5, height = 8)

########################## FIG. 3 - SPECIES DIVERSITY  INTEGRITY ################
######## Fig. 3A beta diverstiy homogeneity #####
# create strip level matrix to capture whole community
head(sp.strip); dim(sp.strip)
# remove habitat info and create species matrix
matrix <-  sp.strip[, -c(1:2,4:6,128)]; rownames(matrix)<- matrix[,1] ; matrix <- matrix[,-c(1)]; matrix

# and point-level
head(sp.point); dim(sp.point)
matrix.point <-  sp.point[, -c(1:3,5:7,129:131)]; rownames(matrix.point)<- matrix.point[,1] ; matrix.point <- matrix.point[,-1]; matrix.point

## basic whittaker index- This is equal to Sorensen index (binary Bray-Curtis in vegan)
bt.div <- betadiver(matrix[colSums(matrix)>0], method = "w");bt.div 

# are there differences between the groups? adonis - yes (use latest version of vegan)
adonis(bt.div ~ summ.strip$habitat.3, method="bray")

# calculate multivariate dispersions (variances; average distance to centroids)
# quickly plot to see dispersions, and differences in composition between groups
bt.disp <- betadisper(bt.div, group = summ.strip$habitat.3, type = "median", bias.adjust = F); plot(bt.disp)
TukeyHSD(bt.disp) # p-value is not significant meaning that group dispersions are homogenous

### principal coordinate analyses https://rstudio-pubs-static.s3.amazonaws.com/78315_a0b12a14d82e47c195b3f2fd506ba6bd.html
# Principal coordinates analysis   pco {ecodist}
library(ecodist)
pcoa.bc <- pco(bt.div) #equivalent to betadisper

# prep df for ggplot 
labs <- paste("PCoA", 1:2, " (", round(100*pcoa.bc$values[1:2]), "%)", sep = ""); labs

#plot.df <-  data.frame(strip.id=rownames(bt.disp$vectors), PCoA1=bt.disp$vectors[,1],PCoA2=bt.disp$vectors[,2] )
plot.df <-  data.frame(strip.id=rownames(bt.disp$vectors), PCoA1=pcoa.bc$vectors[,1],PCoA2=pcoa.bc$vectors[,2] )
plot.df$habitat <- summ.strip$habitat[match(plot.df$strip.id, summ.strip$strip_id)]
plot.df$habitat.3 <-summ.strip$habitat.3[match(plot.df$strip.id, summ.strip$strip_id)]
plot.df$habitat.3 <- factor(plot.df$habitat.3, levels=c("primary", "UPA.Str", "UPA.Int"), labels=c("Primary", "Riparian reserve", "Interior logged"))

fig.3a <- ggscatter( plot.df, x="PCoA1" , y="PCoA2", fill="habitat.3",  shape = "habitat.3",
                      palette = "jco",color="habitat.3",ellipse = T, ellipse.type = "convex",
                      mean.point = T,ellipse.border.remove = T, ellipse.alpha = .25,
                      star.plot = T, size = 3, mean.point.size = 8, star.plot.lwd = .4, star.plot.lty = 3)+
  #stat_ellipse(type = "t",geom="polygon", aes(fill = habitat.3),  alpha = 0.2, show.legend = FALSE, level = 0.75, lty=2, size=.2)+
  scale_fill_manual(name="Forest type",  values=c("#107a7c" , "#56B4E9","#CC79A7"))+
  scale_color_manual(name="Forest type",  values=c("black","black","black"))+
  scale_shape_manual(name="Forest type", values=c(22,21,24))+
  geom_point(data=subset(plot.df, habitat.3=="Primary"), 
             aes( x=PCoA1 , y=PCoA2), color="black", size=6, shape=c(24,21,24,21))+
  ylab(labs[2])+xlab(labs[1])+ theme_bw()+ 
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1.1),
        panel.background = element_blank(), legend.position = "top",axis.text =element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.title=element_blank(), legend.text=element_text(size=16)); fig.3a3

######## Fig. 3B dissimilarity OF RESERVES VS INTERIOR to  primary forest #### 
#### point level ####
#bray.clust <-  hclust(bray.dist, method="complete"); plot(bray.clust)
betadiver(help=TRUE)
bray.dist <-vegdist(matrix.point[,-c(122:124)], method="bray"); bray.dist

# bray disimilarity of pairs of strips, get first pop of the comparison
bray.long <-  gather(as.data.frame(as.matrix(bray.dist)), key=pop1, value=bray);bray.long
# second pop of the comparison in order
bray.long$pop2 <- rep(colnames(as.matrix(bray.dist)),nrow((as.matrix(bray.dist))));bray.long
bray.long$habitat.3.pop2 <- summ.point$habitat.3[match(bray.long$pop2, summ.point$point_id)];bray.long
bray.long$habitat.3.pop1 <- summ.point$habitat.3[match(bray.long$pop1, summ.point$point_id)];bray.long
bray.long$habitat.pop2 <- summ.point$habitat[match(bray.long$pop2, summ.point$point_id)];bray.long
bray.long$habitat.pop1 <- summ.point$habitat[match(bray.long$pop1, summ.point$point_id)];bray.long
# create variable area:site:point (without transect, to draw lines between points), e.g. UPA2p_3
bray.long$area.site.point.pop1 <- paste(substr(bray.long$pop1,0,6), substr(bray.long$pop1,8,8), sep="") ;bray.long
bray.long$area.site.point.pop2 <- paste(substr(bray.long$pop2,0,6), substr(bray.long$pop2,8,8), sep="") ;bray.long
bray.long$area.site.pop1 <- paste(substr(bray.long$pop1,0,5), sep="") ;bray.long
bray.long$area.site.pop2 <- paste(substr(bray.long$pop2,0,5), sep="") ;bray.long

# add what they were compared to
bray.long$area.site.point.pop1.area.site.point.pop2 <- paste(bray.long$area.site.point.pop1,bray.long$area.site.point.pop2, sep="_") ;bray.long

# subset so that one of the pops is always primary (either interior or stream)
bray.long.prim <- subset(bray.long, habitat.3.pop2=="primary");bray.long.prim
head(bray.long.prim)

# are comparisons prim of the same area.sites ('adjacent')? i.e. PRIMd vs PRIMd
bray.long.prim$prim.vs.prim.same.site <-ifelse(bray.long.prim$area.site.pop1==bray.long.prim$area.site.pop2, "yes", "no"); head(bray.long.prim)

# summarise mean difference of each site to all primary sites (stream and interior separately)
bray.long.prim <- subset(bray.long.prim,bray!=0); head(bray.long.prim)
bray.long.prim.mean <-  summarise(group_by(bray.long.prim, pop1, habitat.pop1, habitat.pop2, habitat.3.pop1, 
                                           habitat.3.pop2, area.site.point.pop1, prim.vs.prim.same.site),
                                  mean.bray.to.prim=mean(bray)); bray.long.prim.mean
# means for text
summarise(group_by(subset(bray.long.prim.mean, habitat.3.pop1!=""&prim.vs.prim.same.site=="no"), habitat.3.pop1),
          n=n(), mean.bray.to.prim.mean=mean(mean.bray.to.prim), sd.bray.to.prim=sd(mean.bray.to.prim))

# primary sites non adjacent
bray.long.prim.mean$habitat.3.pop1 <- factor(bray.long.prim.mean$habitat.3.pop1, levels=c("primary", "UPA.Str", "UPA.Int"), labels=c("Primary\nnon-adjacent", "Riparian\nreserve", "Interior\nlogged"))
bray.long.prim$habitat.3.pop1 <- factor(bray.long.prim$habitat.3.pop1, levels=c("primary", "UPA.Str", "UPA.Int"), labels=c("Primary", "Riparian\nreserve", "Interior\nlogged"))
bray.long.prim.mean$site.pop1 <- substr(bray.long.prim.mean$pop1,0,7); bray.long.prim.mean$site.pop1


# add colours to lines between rip and int, does dissimilarity increase between pairs of closest points?  
bray.long.prim.mean$area.site.point.pop1_habitat.pop2 <- paste(bray.long.prim.mean$area.site.point.pop1, bray.long.prim.mean$habitat.pop2, sep = "_")
change.diss <- subset(bray.long.prim.mean, habitat.3.pop1!=""&prim.vs.prim.same.site=="no")[,c("habitat.pop1", "area.site.point.pop1_habitat.pop2",
                                                                                        "mean.bray.to.prim")] %>%
  tidyr::spread(habitat.pop1 , mean.bray.to.prim); change.diss <- change.diss[,-c(2,3)]; head(change.diss) 

change.diss$increase.from.upa.str_upa.int <- if_else(change.diss$UPA.Int>change.diss$UPA.Str, "yes", "no"); head(change.diss) 
bray.long.prim.mean$increase.from.upa.str_upa.int <- change.diss$increase.from.upa.str_upa.int[match(bray.long.prim.mean$area.site.point.pop1_habitat.pop2,
                                                                                                     change.diss$area.site.point.pop1_habitat.pop2)]; head(bray.long.prim.mean)
change.diss$increase.val.from.upa.str_upa.int <- change.diss$UPA.Int-change.diss$UPA.Str; head(change.diss) 
bray.long.prim.mean$increase.val.from.upa.str_upa.int <- change.diss$increase.val.from.upa.str_upa.int[match(bray.long.prim.mean$area.site.point.pop1_habitat.pop2,
                                                                                                     change.diss$area.site.point.pop1_habitat.pop2)]; head(bray.long.prim.mean)

# means for text
mean(subset(bray.long.prim.mean, habitat.3.pop1!=""&prim.vs.prim.same.site=="no"&increase.from.upa.str_upa.int!=""&increase.from.upa.str_upa.int=="yes")$increase.val.from.upa.str_upa.int)
length(subset(bray.long.prim.mean, habitat.3.pop1!=""&prim.vs.prim.same.site=="no"&increase.from.upa.str_upa.int!=""&increase.from.upa.str_upa.int=="yes")$increase.val.from.upa.str_upa.int)
mean(subset(bray.long.prim.mean, habitat.3.pop1!=""&prim.vs.prim.same.site=="no"&increase.from.upa.str_upa.int!=""&increase.from.upa.str_upa.int=="no")$increase.val.from.upa.str_upa.int)
length(subset(bray.long.prim.mean, habitat.3.pop1!=""&prim.vs.prim.same.site=="no"&increase.from.upa.str_upa.int!=""&increase.from.upa.str_upa.int=="no")$increase.val.from.upa.str_upa.int)

# plot
my_comparisons <- list(c("Primary\nnon-adjacent", "Riparian\nreserve"),c("Interior\nlogged", "Riparian\nreserve"), c("Primary\nnon-adjacent", "Interior\nlogged") )
my_comparisons1 <- list(c("Interior\nlogged", "Riparian\nreserve") )

fig3.b <- ggplot(aes(x=habitat.3.pop1 , y=mean.bray.to.prim, fill=habitat.3.pop1 ), 
                 data=subset(bray.long.prim.mean, habitat.3.pop1!=""&prim.vs.prim.same.site=="no"))+
  geom_boxplot(outlier.shape = NA)+
  annotate(geom = "rect",xmin = 2.01, xmax = 2.99, ymin = -Inf, ymax=Inf, color="white", fill="white")+
  scale_fill_manual(name="Forest type",  values=c("#107a7c" , "#56B4E9","#CC79A7"))+
  scale_color_manual(name="Forest type",  values=c("#107a7c" ,"#56B4E9","#CC79A7"))+
  #stat_compare_means(method = "t.test", comparisons = my_comparisons,label = "p.signif", paired=F, tip.length = 0.0,vjust = .1 , step.increase = 0.05, size=5) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons1,label = "p.signif", paired=T,
                     tip.length = 0.0,vjust = .5 , step.increase = 0.01, size=8) +
  #scale_shape_manual(name="Forest type", values=c(22,21,24))+
  geom_jitter(aes(group=habitat.pop2, shape=habitat.pop2),width=.05, alpha=.9,size=3, colour="black", fill=NA)+
  scale_shape_manual(values=c(21,24))+
  geom_line(inherit.aes=F, data=subset(bray.long.prim.mean, habitat.3.pop1!=""&prim.vs.prim.same.site=="no"&increase.from.upa.str_upa.int!=""), 
            aes(x=habitat.3.pop1 , y=mean.bray.to.prim, group=area.site.point.pop1_habitat.pop2, colour=increase.from.upa.str_upa.int), alpha=.7)+
  scale_color_manual(name="increase",  values=c("darkgrey","orange"))+
  annotate(geom = "text",x = 1, y = 0.535, label="mean change= +0.05 (n=54)", color="orange", size=7,hjust = 0)+
  annotate(geom = "text",x = 1, y = 0.52, label="mean change= -0.02 (n=42)", color="grey", size=7,hjust = 0)+
  ylim(0.52,0.81)+
  #EnvStats::stat_n_text()+
  ylab("Mean dissimilarity to primary points (Bray-Curtis)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 16),
        axis.text =element_text(size = 16), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_text(size=17, face = c("plain", "bold","plain")),
        axis.ticks.x = element_line(colour = c("black","blue","black"), size = c(1,7,1))); fig3.b

#### t-test and means ####
t.test(data=subset(bray.long.prim.mean, habitat.3.pop1!=""&habitat.3.pop1!="Primary\nnon-adjacent"),  
       mean.bray.to.prim~ habitat.3.pop1, 
       paired = T)
summarise(group_by(subset(bray.long.prim.mean, habitat.3.pop1!=""&habitat.3.pop1!="Primary\nnon-adjacent"), habitat.3.pop1),
          mean=mean(mean.bray.to.prim),
          n=n())

######## Fig. 3C analyses - disimilarities against canopy height #####
## check if any environmental variables explain "primaryness"
## canopy height makes sites more "primary" like
bray.long.prim.mean$max.canopy <- veg.point$max.canopy[match(bray.long.prim.mean$pop1, veg.point$point_code)]
bray.long.prim.mean$canopy.cover<- veg.point$canopy.cover[match(bray.long.prim.mean$pop1, veg.point$point_code)]
bray.long.prim.mean$mean.canopy<- veg.point$mean.canopy[match(bray.long.prim.mean$pop1, veg.point$point_code)]
bray.long.prim.mean$diam.tree.total<- veg.point$diam.tree.total[match(bray.long.prim.mean$pop1, veg.point$point_code)]
bray.long.prim.mean$diam.tree.mean<- veg.point$diam.tree.mean[match(bray.long.prim.mean$pop1, veg.point$point_code)]
bray.long.prim.mean$stream.width<- veg.point$stream.width[match(bray.long.prim.mean$pop1, veg.point$point_code)]
bray.long.prim.mean$upa <- substr(bray.long.prim.mean$pop1, 0,4)

## model if environmental variables change dissimilarity ###
## primary
head(bray.long.prim.mean)
summary(lm(mean.bray.to.prim~canopy.cover+diam.tree.total+stream.width+mean.canopy, 
           data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Str")))
lm1 <- lm(mean.bray.to.prim~canopy.cover+diam.tree.total+stream.width+mean.canopy,
          data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Str")); summary(lm1)
library(relaimpo); calc.relimp(lm1)
lm1 <- lm(mean.bray.to.prim~mean.canopy,
          data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Str")); summary(lm1)


# use beta regression
library(betareg)
br1 <- betareg(mean.bray.to.prim ~canopy.cover+diam.tree.mean*mean.canopy +stream.width,
          data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Str")); summary(br1)
brtest::lrtest(br1); plot(br1)
br1 <- betareg(mean.bray.to.prim~mean.canopy,
          data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Str")); summary(br1)
lmtest::lrtest(br1)

## UPA.int
br2 <-betareg(mean.bray.to.prim~canopy.cover+diam.tree.mean*mean.canopy, 
       data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Int")); summary(br2)
lmtest::lrtest(br2)

#### glmer dissimilarity explained by E vairables ####
bray.long.prim.mean$area2 <- if_else(substr(bray.long.prim.mean$area.site.point.pop1,0,1)=="P",
                                     substr(bray.long.prim.mean$area.site.point.pop1,0,5),
                                     substr(bray.long.prim.mean$area.site.point.pop1,0,4) ); unique(bray.long.prim.mean$area2)

#install.packages("glmmTMB", type="source")
library(glmmTMB)
library(DHARMa)
library(effects)
library(ggplot2)
library(broom.mixed)

# scale variables
bray.long.prim.mean$mean.canopy.scaled <- scale(bray.long.prim.mean$mean.canopy, center = T)
bray.long.prim.mean$canopy.cover.scaled <- scale(bray.long.prim.mean$canopy.cover, center = T)
bray.long.prim.mean$diam.tree.mean.scaled <- scale(bray.long.prim.mean$diam.tree.mean, center = T)
bray.long.prim.mean$canopy.cover.scaled <- scale(bray.long.prim.mean$canopy.cover, center = T)
bray.long.prim.mean$stream.width.scaled <- scale(bray.long.prim.mean$stream.width, center = T)

#### riparian reserves points dissimilarity #####
# use site (n=6) to control for spatial autocorrelation
# mean tree diameter and mean canopy height interact
glmmTMB.str <- glmmTMB(mean.bray.to.prim~ canopy.cover.scaled + stream.width.scaled +
                       diam.tree.mean.scaled *mean.canopy.scaled + (1 | site.pop1), 
                       data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Str"&mean.canopy!=""), 
                     family="beta_family"); summary(glmmTMB.str)

ae <- allEffects(glmmTMB.str); plot(ae)
car::Anova(glmmTMB.str,type="III")

res = simulateResiduals(glmmTMB.str); plot(res, rank = T)
options(scipen = 999)
broom.mixed::tidy(glmmTMB.str)

#### logged interior points dissimilarity #####
library(car)
glmmTMB.int <- glmmTMB(mean.bray.to.prim~ canopy.cover.scaled  + 
                          diam.tree.mean.scaled *mean.canopy.scaled+ (1 | site.pop1),
                        data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Int"&mean.canopy!=""),
                        family="beta_family", na.action = "na.fail"); summary(glmmTMB.int)

ae <- allEffects(glmmTMB.int); plot(ae)
car::Anova(glmmTMB.int,type="III")

res = simulateResiduals(glmmTMB.int); plot(res, rank = T)
broom.mixed::tidy(glmmTMB.int)

## drop nothing
glmmTMB.int_d1 <- drop1(glmmTMB.int,test="Chisq"); print(glmmTMB.int_d1)
aa <- augment(glmmTMB.int, data=subset(bray.long.prim.mean, habitat.pop1=="UPA.Int"&mean.canopy!=""))
ggplot(aa, aes(.fitted,.resid)) + geom_line(aes(group=site.pop1),colour="gray") + geom_point(aes(colour=upa))   + geom_smooth()



######## Fig. 3C. plot - get primary interior mean bray and mean canopy, with sd #####
prim.df <-  data.frame(bray.mean=mean(subset(bray.long.prim.mean, habitat.pop1=="PRI.Int")$mean.bray.to.prim),
                       bray.sd=sd(subset(bray.long.prim.mean, habitat.pop1=="PRI.Int")$mean.bray.to.prim),
                       mean.canopy.mean=mean(subset(bray.long.prim.mean, habitat.pop1=="PRI.Int")$mean.canopy),
                       mean.canopy.sd=sd(subset(bray.long.prim.mean, habitat.pop1=="PRI.Int")$mean.canopy)); prim.df

bray.long.prim.mean$habitat.pop2 <- factor(bray.long.prim.mean$habitat.pop2, levels=c("PRI.Str", "PRI.Int"), labels=c("Primary\nriparian", "Primary\ninterior"))

fig3.c <-  ggplot(aes(x=mean.canopy  , y=mean.bray.to.prim,  fill=habitat.pop1,  shape = habitat.pop1), data=subset(bray.long.prim.mean, habitat.3.pop1!="Primary\nnon-adjacent"))+
  facet_wrap(~habitat.3.pop1, nrow = 2)+ 
  geom_vline(xintercept = prim.df$mean.canopy.mean, lty=2, color = "#107a7c" )+
  stat_smooth(aes(x=mean.canopy  , y=mean.bray.to.prim,fill=habitat.pop1,colour=habitat.pop1), data=subset(bray.long.prim.mean, habitat.3.pop1=="Riparian\nreserve"),method="lm", inherit.aes = F)+
  #geom_line(aes(y = predict(gy_logit, GasolineYield),   colour = "logit", linetype = "logit")) +
  stat_cor(aes(colour=habitat.pop1, label =paste(..rr.label.., cut(..p.., breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
              labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),label.x.npc = 0.45, label.y = 0.75, inherit.aes = T, size=6)+ 
  geom_point(size=3)+
  scale_fill_manual(name="Forest type",  values=c("#CC79A7","#56B4E9"))+
  scale_color_manual(name="Forest type",  values=c("#CC79A7","#56B4E9"))+
  scale_shape_manual(name="Forest type", values=c(24,21))+
  scale_x_continuous(limits=c(8,25) , expand = c(0, 0)) +
  scale_y_continuous(limits=c(0.54,0.77), expand = c(0, 0)) +
  ylab("Mean dissimilarity to primary points (Bray-Curtis)")+xlab("Canopy height (m)")+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=18),
        axis.text =element_text(size = 16), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()); fig3.c

# divide by whether is comparison to primary interior points or primary riparian points
fig3.c.SI <-  ggplot(aes(x=mean.canopy  , y=mean.bray.to.prim,  fill=habitat.pop1,  shape = habitat.pop1), data=subset(bray.long.prim.mean, habitat.3.pop1!="Primary\nnon-adjacent"))+
  facet_grid(habitat.3.pop1~habitat.pop2)+ 
  geom_vline(xintercept = prim.df$mean.canopy.mean, lty=2, color = "#107a7c" )+
  stat_smooth(aes(x=mean.canopy  , y=mean.bray.to.prim,fill=habitat.pop1,colour=habitat.pop1), data=subset(bray.long.prim.mean, habitat.3.pop1=="Riparian\nreserve"),method="lm", inherit.aes = F)+
  #geom_line(aes(y = predict(gy_logit, GasolineYield),   colour = "logit", linetype = "logit")) +
  stat_cor(aes(colour=habitat.pop1, label =paste(..rr.label.., cut(..p.., breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                   labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),label.x.npc = 0.45, label.y = 0.75, inherit.aes = T, size=6)+ 
  geom_point(size=3)+
  scale_fill_manual(name="Forest type",  values=c("#CC79A7","#56B4E9"))+
  scale_color_manual(name="Forest type",  values=c("#CC79A7","#56B4E9"))+
  scale_shape_manual(name="Forest type", values=c(24,21))+
  scale_x_continuous(limits=c(8,25) , expand = c(0, 0)) +
  scale_y_continuous(limits=c(0.54,0.77), expand = c(0, 0)) +
  ylab("Mean dissimilarity to primary points (Bray-Curtis)")+xlab("Canopy height (m)")+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text =element_text(size = 16), 
        strip.background = element_blank(),
        strip.text.x  = element_text(size = 16),        strip.text.y  = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()); fig3.c.SI
ggsave("plots/SI/fig3.SI.png", width = 8, height = 8)

#### compose figure 3 ####
cowplot::plot_grid(fig.3a3, fig3.b,  fig3.c, ncol = 3,align = "hv", rel_widths = c(.9,.8,.6), 
                   labels = c("A","B", "C"), scale = c(.9,.95,.95), label_size = 24)

#ggsave("plots/fig3.png", width = 16, height = 7)


########################## FIG. 4 PHYLOGENY ##########################
# make matrix into long format
sp.habitat.long <- gather(sp.habitat,genus_sp,abundance,-location,-habitat); head(sp.habitat.long ) ;sp.habitat.long <- sp.habitat.long[,-1]

# wide so that habitats are columns
sp.habitat.long.w <- tidyr::spread(sp.habitat.long,habitat,abundance, fill = F)

# create one primary with both interior and riparian
sp.habitat.long.w$primary <-sp.habitat.long.w$PRI.Int+sp.habitat.long.w$PRI.Str

# normalise by number of points in each habitat type; (i.e. per point(site), per sp)
summarise(group_by(sp, habitat), n.points=length(unique(point_id)))
sp.habitat.long.w$primary.norm <- sp.habitat.long.w$primary/16
sp.habitat.long.w$UPA.Int.norm <- sp.habitat.long.w$UPA.Int/24
sp.habitat.long.w$UPA.Str.norm <- sp.habitat.long.w$UPA.Str/24; head(sp.habitat.long.w)
sp.habitat.norm <- sp.habitat.long.w[,c(1,7,8,9)]
sp.habitat.norm$genus_sp <- sp$genus_sp[match(sp.habitat.norm$genus_sp, sp$genus_sp)]; 
head(sp.habitat.norm)

# how many primary species were absent from both reserves and logged forests? /120
head(sp.habitat.long.w)
nrow(subset(sp.habitat.long.w, primary.norm!=0 & UPA.Int.norm==0 & UPA.Str.norm==0 ))


# how are butterflies distributed across the three habitats? percetanges of total abundance
sp.habitat.norm$total.ab.norm <- sp.habitat.norm$primary.norm+sp.habitat.norm$UPA.Int.norm+sp.habitat.norm$UPA.Str.norm
sp.habitat.norm$primary.norm.perc <- (sp.habitat.norm$primary.norm/sp.habitat.norm$total.ab.norm)*100
sp.habitat.norm$UPA.Int.norm.perc <- (sp.habitat.norm$UPA.Int.norm/sp.habitat.norm$total.ab.norm)*100
sp.habitat.norm$UPA.Str.norm.perc <- (sp.habitat.norm$UPA.Str.norm/sp.habitat.norm$total.ab.norm)*100; head(sp.habitat.norm)

# normalised % abundance in riparian/logged - ab in primary (negative means loss compared to primary)
sp.habitat.norm$upa.int.minus.prim <- round(sp.habitat.norm$UPA.Int.norm.perc-sp.habitat.norm$primary.norm.perc)
sp.habitat.norm$upa.str.minus.prim <- round(sp.habitat.norm$UPA.Str.norm.perc-sp.habitat.norm$primary.norm.perc)
head(sp.habitat.norm)
sp.habitat.norm$ab.loss.in.logged <- if_else(sp.habitat.norm$upa.int.minus.prim==-100, "lost",  
                                             if_else(sp.habitat.norm$upa.int.minus.prim==0,"equal", 
                                                     if_else(sp.habitat.norm$upa.int.minus.prim<0,"loss.ab", "gain.ab" )) )
sp.habitat.norm$ab.loss.in.reserve <- if_else(sp.habitat.norm$upa.str.minus.prim==-100, "lost",  
                                             if_else(sp.habitat.norm$upa.str.minus.prim==0,"equal", 
                                                     if_else(sp.habitat.norm$upa.str.minus.prim<0,"loss.ab", "gain.ab" )) )



### plot data on tree ####
# grab diff in perc ab
head(sp.habitat.norm)
sp.habitat.df.loss <-  sp.habitat.norm[,c("genus_sp", "upa.int.minus.prim","upa.str.minus.prim")]; head(sp.habitat.df.loss)

# add the obvious, no difference between primary and primary
sp.habitat.df.loss$prim <- rep(0, nrow(sp.habitat.df.loss))

# make into matrix for ggtree to read
row.names(sp.habitat.df.loss)<- as.character(sp.habitat.df.loss$genus_sp)
sp.habitat.df.loss <- sp.habitat.df.loss[,c(4,3,2)]; head(sp.habitat.df.loss)
names(sp.habitat.df.loss)<- c("Primary", "Reserve", "Logged")
head(sp.habitat.df.loss); str(sp.habitat.df.loss)

# change species names (remove _) so that they match pretty tree
row.names(sp.habitat.df.loss)<- gsub("_", " ", row.names(sp.habitat.df.loss))

# plot
p <- gheatmap(p1, abs(sp.habitat.df.loss), colnames=T,offset=10.5, width=0.6,color="black",
              colnames_angle=45,colnames_offset_y = 0.5,colnames_offset_x = .8, colnames_position="top")+
  scale_fill_viridis_c( name="abundance loss \n compared to primary")+
  scale_fill_gradient(low = "#107a7c", high = "white",limits=c(0,100),
                      name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 1.3)); p

## add abundance in primary forest as circles
col <- data.frame(genus_sp=tree.genus.sp.prim$tip.label, ab.prim=NA)
sp.habitat.long.w$genus_sp.space <- gsub("_", " ", sp.habitat.long.w$genus_sp )
col$ab.prim <- as.numeric(sp.habitat.long.w$primary[match(col$genus_sp, sp.habitat.long.w$genus_sp.space)])
str(col)
p <- p %<+% col; p + geom_tippoint(aes(x=x-0.01, size=ab.prim),colour="darkgrey", fill="#107a7c", shape=21) +
  #scale_color_manual(values=c("#107a7c"), guide='none') +
  labs(fill="Abundance primary")+
  ggplot2::xlim(0, 24) 


ggsave("plots/fig4.2.png", height = 12, width = 8)





###### Sophie's tree- make tree of genus/sp not in mena's tree with sophie's tree ######
library(ggtree)
tree.sophie.genus.sp <- read.nexus("1.data/tree.3genes.genus.sp.nex")
plot(tree.sophie.genus.sp)
### prep treeto only show species found at least once in primary forest  #####
prim.sp <- unique(subset(sp.habitat.norm, genus_sp %in% subset(sp.habitat.long.w, primary>0)$genus_sp)$genus_sp)

# prim sp total - 72
length(subset(sp.habitat.long.w, primary>0)$genus_sp)
prim.sp <- unique(subset(sp.habitat.norm, genus_sp %in% subset(sp.habitat.long.w, primary>0)$genus_sp)$genus_sp)

# drop tips that are not prim.sp
tree.sophie.genus.sp.prim <- drop.tip(tree.sophie.genus.sp ,subset(tree.sophie.genus.sp$tip.label, !(tree.sophie.genus.sp$tip.label%in%prim.sp)) ); plot(tree.sophie.genus.sp.prim)

# make nice species names
tree.sophie.genus.sp.prim$tip.label <- gsub("_", " ", tree.sophie.genus.sp.prim$tip.label)
length(tree.sophie.genus.sp.prim$tip.label)

## make tree with all species found in jamari
tree.sophie.genus.sp.all <- drop.tip(tree.sophie.genus.sp ,subset(tree.sophie.genus.sp$tip.label, !(tree.sophie.genus.sp$tip.label%in%sp.habitat.long.w$genus_sp)) ); plot(tree.sophie.genus.sp.all)

######## Fig. 4B plot base tree ####
# could remove branch lenghts  branch.length="none"
p1 <- ggtree(tree.sophie.genus.sp.prim, color="grey", size=1,branch.length="none" ) + 
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3, colour="grey", size=0)  +
  geom_tiplab(align=TRUE,  linetype='dashed', linesize=0,colour="black", size=4.5, fontface='italic', offset = 0.4)+xlim(0,40)  ; p1

### plot data on tree ####
# grab diff in perc ab
head(sp.habitat.norm)
sp.habitat.df.loss <-  sp.habitat.norm[,c("genus_sp", "upa.int.minus.prim","upa.str.minus.prim")]; head(sp.habitat.df.loss)

# add the obvious, no difference between primary and primary
sp.habitat.df.loss$prim <- rep(0, nrow(sp.habitat.df.loss))

# make into matrix for ggtree to read
row.names(sp.habitat.df.loss)<- as.character(sp.habitat.df.loss$genus_sp)
sp.habitat.df.loss <- sp.habitat.df.loss[,c(4,3,2)]; head(sp.habitat.df.loss)
names(sp.habitat.df.loss)<- c("Primary", "Reserve", "Logged")
head(sp.habitat.df.loss); str(sp.habitat.df.loss)

# change species names (remove _) so that they match pretty tree
row.names(sp.habitat.df.loss)<- gsub("_", " ", row.names(sp.habitat.df.loss))

# plot
p <- gheatmap(p1, abs(sp.habitat.df.loss[,2:3]), colnames=T,offset=15.5, width=0.2,color="black",
              colnames_angle=45,colnames_offset_y = 1.5,colnames_offset_x = .8, colnames_position="top")+
  scale_fill_viridis_c( name="abundance loss \n compared to primary")+
  scale_fill_gradient(low = "#107a7c", high = "white",limits=c(0,100),
                      name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 1.3)); p

## add abundance in primary forest as circles
col <- data.frame(genus_sp=tree.sophie.genus.sp.prim$tip.label, ab.prim=NA)
sp.habitat.long.w$genus_sp.space <- gsub("_", " ", sp.habitat.long.w$genus_sp )
col$ab.prim <- as.numeric(sp.habitat.long.w$primary[match(col$genus_sp, sp.habitat.long.w$genus_sp.space)])
str(col)
p <- p %<+% col; p + geom_tippoint(aes(x=x-0.01, size=ab.prim),colour="darkgrey", fill="#107a7c", shape=21) +
  #scale_color_manual(values=c("#107a7c"), guide='none') +
  labs(fill="Abundance primary")+
  ggplot2::xlim(0, 34) 

#ggsave("plots/fig4.png", height = 12, width = 8)

## summarise by genus and plot change with habitat ###

#### normalised abundances by point , for plotting per species ####
head(sp.habitat.long.w) 
sp.habitat.norm.plot <-  sp.habitat.norm[,c("genus_sp","primary.norm.perc","UPA.Int.norm.perc","UPA.Str.norm.perc","upa.int.minus.prim", "upa.str.minus.prim")]; head(sp.habitat.norm.plot)
sp.habitat.norm.plot$prim.minus.prim <- 0
sp.habitat.norm.plot.long <- gather(sp.habitat.norm.plot[,c("genus_sp","primary.norm.perc","UPA.Int.norm.perc","UPA.Str.norm.perc")], habitat, norm.ab.perc, -genus_sp); head(sp.habitat.norm.plot.long)
sp.habitat.norm.plot.long$diff.to.prim <-NA
sp.habitat.norm.plot.long[sp.habitat.norm.plot.long$habitat=="primary.norm.perc",]$diff.to.prim <- sp.habitat.norm.plot$prim.minus.prim[match(sp.habitat.norm.plot$genus_sp,sp.habitat.norm.plot$genus_sp)]; head(sp.habitat.norm.plot.long)
sp.habitat.norm.plot.long[sp.habitat.norm.plot.long$habitat=="UPA.Int.norm.perc",]$diff.to.prim <- sp.habitat.norm.plot$upa.int.minus.prim[match(sp.habitat.norm.plot$genus_sp,sp.habitat.norm.plot$genus_sp)]; head(sp.habitat.norm.plot.long)
sp.habitat.norm.plot.long[sp.habitat.norm.plot.long$habitat=="UPA.Str.norm.perc",]$diff.to.prim <- sp.habitat.norm.plot$upa.str.minus.prim[match(sp.habitat.norm.plot$genus_sp,sp.habitat.norm.plot$genus_sp)]; head(sp.habitat.norm.plot.long)

sp.habitat.norm.plot.long$habitat <-factor(sp.habitat.norm.plot.long$habitat, 
                                        levels=c("primary.norm.perc", "UPA.Str.norm.perc","UPA.Int.norm.perc"),
                                        labels=c("Primary", "Riparian\nreserve", "Interior\nlogged"))

sp.habitat.long.w$total.abundance <- sp.habitat.long.w$UPA.Int +sp.habitat.long.w$UPA.Str +sp.habitat.long.w$primary
sp.habitat.norm.plot.long$genus <- unlist(lapply(str_split(sp.habitat.norm.plot.long$genus_sp, "_"), head, 1))
  
######## Fig 4A. plot change from rip to log ######
sp.habitat.norm.plot.long$diff.to.prim.abs <- abs(sp.habitat.norm.plot.long$diff.to.prim)
#&diff.to.prim.abs<100
plot.change.df <- subset(sp.habitat.norm.plot.long ,habitat!= "Primary"&genus_sp%in%prim.sp)[,c("genus_sp", "habitat", "diff.to.prim.abs")] %>%
  tidyr::spread(habitat, diff.to.prim.abs) %>%
  `colnames<-`(c("genus_sp", "rip", "log")) %>%
  subset((rip!=100&log!=100))%>% #exclude those that disappear
  dplyr::mutate(is_increasing = rip < log) %>%
  tidyr::gather("habitat", "diff.to.prim.abs", 2:3); plot.change.df

plot.change.df$habitat <-factor(plot.change.df$habitat,  levels=c("rip", "log"), labels=c("Riparian\nreserve", "Interior\nlogged"))
plot.change.df$subfamily <- sp$subfamily[match(plot.change.df$genus_sp, sp$genus_sp)]
plot.change.df$tribe <- sp$tribe[match(plot.change.df$genus_sp, sp$genus_sp)]

compare_means(data=subset(plot.change.df),  diff.to.prim.abs~ habitat, 
              paired = TRUE, method="t.test", label = "p.format", size=6)

t.test(data=subset(plot.change.df),  diff.to.prim.abs~ habitat, 
       paired = TRUE)

## get means and numbers of species 
plot.change.df.wide <- spread(plot.change.df,  habitat, diff.to.prim.abs); head(plot.change.df.wide)
plot.change.df.wide$diff.int.log.minus.rip.res <-  plot.change.df.wide$`Interior\nlogged` -plot.change.df.wide$`Riparian\nreserve`
plot.change.df.wide$habitat.more.similar.to.primary <- if_else(plot.change.df.wide$is_increasing==TRUE, "Riparian\nreserve", "Interior\nlogged"); head(plot.change.df.wide)

mean(subset(plot.change.df.wide, is_increasing ==TRUE )$diff.int.log.minus.rip.res)
length(subset(plot.change.df.wide, is_increasing ==TRUE )$diff.int.log.minus.rip.res)

mean(subset(plot.change.df.wide, is_increasing ==FALSE )$diff.int.log.minus.rip.res)
length(subset(plot.change.df.wide, is_increasing ==FALSE  )$diff.int.log.minus.rip.res)

### plot ####
# add mean changes and number of species 
my_comparisons1 <- list(c("Interior\nlogged", "Riparian\nreserve") )

ggplot(data=subset(plot.change.df), aes(x = habitat, y = diff.to.prim.abs)) +
  geom_boxplot(fill=c("#56B4E9", "#CC79A7"), alpha=0.6)+
  ylab("% abundance change with respect to primary, per species")+
  annotate(geom = "rect",xmin = 1.01, xmax = 1.99, ymin = -Inf, ymax=Inf, color="white", fill="white")+
  geom_point(size=5, shape=21, colour="black", alpha=0.8, aes(fill=diff.to.prim.abs, colour=diff.to.prim.abs)) +
  scale_fill_gradient(low = "#107a7c", high = "white", limits=c(0,100), 
                      name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  scale_colour_gradient(low = "#107a7c", high = "white", limits=c(0,100), 
                        name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  geom_line(aes(group = genus_sp, col = is_increasing), size=.5, alpha=0.8) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons1,label = "p.signif", paired=T,
                     tip.length = 0.0,vjust = .5 , step.increase = 0.01, size=8) +
  # annotate(geom = "text",x = 0.5, y = -5, label="Int. log. > ab. change, mean= +27 (n=36)", color="orange", size=6,hjust = 0)+
  # annotate(geom = "text",x = 0.5, y = -10, label="Int. log. < ab. change, mean= -14 (n=29)", color="grey40", size=6,hjust = 0)+
  annotate(geom = "text",x = 0.5, y = -5, label="mean change= +27 (n=36)", color="orange", size=7,hjust = 0)+
  annotate(geom = "text",x = 0.5, y = -10, label="mean change= -14 (n=29)", color="grey40", size=7,hjust = 0)+
  scale_colour_manual(values = c( "grey40", "orange"))+
  ylim(c(-10,95))+
  #stat_compare_means(paired = TRUE, method="t.test", label = "p.format", size=6)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 18),
        axis.text.y =element_text(size = 16), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=20, face = c("bold","plain")),
        axis.ticks.x = element_line(colour = c("blue","black"), size = c(7,1)))
  
ggsave("plots/fig4/4.b.change.rip.log.prim.sp.png", height = 8, width = 5.5 )

### plot fig 4 C new one for resubmission, too complicated
# spread difference by habitat
plot.change.df.wide <- spread(plot.change.df,  habitat, diff.to.prim.abs); head(plot.change.df.wide)

plot.change.df.wide$diff.int.log.minus.rip.res <-  plot.change.df.wide$`Interior\nlogged` -plot.change.df.wide$`Riparian\nreserve`
plot.change.df.wide$habitat.more.similar.to.primary <- if_else(plot.change.df.wide$is_increasing==TRUE, "Riparian\nreserve", "Interior\nlogged"); head(plot.change.df.wide)


mean(subset(plot.change.df.wide, is_increasing ==TRUE )$`Riparian\nreserve`)
length(subset(plot.change.df.wide, is_increasing ==TRUE )$diff.int.log.minus.rip.res)

mean(subset(plot.change.df.wide, is_increasing ==FALSE )$diff.int.log.minus.rip.res)
length(subset(plot.change.df.wide, is_increasing ==FALSE  )$diff.int.log.minus.rip.res)

plot.change.df.wide.summ <- summarise(group_by(plot.change.df.wide, is_increasing,habitat.more.similar.to.primary),
                                      n=n(),
                                      diff.to.prim.in.best=if_else(habitat.more.similar.to.primary=="Interior\nlogged", mean(`Interior\nlogged`), mean(`Riparian\nreserve`)),
                                      diff.int.log.minus.rip.res.mean=mean(diff.int.log.minus.rip.res)); plot.change.df.wide.summ 

plot.change.df.wide.summ <-plot.change.df.wide.summ[!duplicated(plot.change.df.wide.summ),]


# how much better off are they depending on whether a species does better in riparian vs logged?
plot.change.df.wide.summ$habitat.more.similar.to.primary <-factor(plot.change.df.wide.summ$habitat.more.similar.to.primary ,  level=c("Riparian\nreserve", "Interior\nlogged")); plot.change.df.wide.summ$habitat.more.similar.to.primary 
library(EnvStats)

ggplot(data=subset(plot.change.df.wide.summ ), aes(x = factor(habitat.more.similar.to.primary, level=c("Riparian\nreserve", "Interior\nlogged"),
                                                              labels=c("Riparian reserves more\nsimilar to primary", "Logged interiors more\nsimilar to primary")), y = abs(diff.int.log.minus.rip.res.mean), fill=habitat.more.similar.to.primary)) +
  geom_point( data=plot.change.df.wide, aes(y=abs(diff.int.log.minus.rip.res), x=factor(habitat.more.similar.to.primary, level=c("Riparian\nreserve", "Interior\nlogged"),
                                                                                        labels=c("Riparian reserves more\nsimilar to primary", "Logged interiors more\nsimilar to primary"))), shape=21) +
  geom_point(aes(size=n), alpha=0.6, shape=21)+
  scale_size_continuous(range=c(20,40))+
  ylab("Difference in abundance between\nriparian reserves and interior logged")+
  xlab("Species with abundance more similar to primary forest")+
  scale_fill_manual(values = c( "grey40", "orange"))+
  stat_n_text(y.pos = -5,size=6, inherit.aes =F, data=plot.change.df.wide, aes(y=abs(diff.int.log.minus.rip.res), x=factor(habitat.more.similar.to.primary, level=c("Riparian\nreserve", "Interior\nlogged"),
                                                                                                                           labels=c("Riparian reserves more\nsimilar to primary", "Logged interiors more\nsimilar to primary"))))+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 16),
        axis.text =element_text(size = 16), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=16, face = c("bold","plain")),
        axis.ticks.x = element_line(colour = c("blue","black"), size = c(7,1)))


######## Fig 4C plot several species per genus ######
genus.n <- summarise(group_by(sp, genus),
          n=n()); genus.n


taygetis.p <- ggplot(subset(sp.habitat.norm.plot.long , genus=="Taygetis"&genus_sp%in%prim.sp), 
                 aes(x=habitat, y=norm.ab.perc, fill=abs(diff.to.prim), colour=abs(diff.to.prim)))+
  geom_line(aes(x=habitat, y=norm.ab.perc,group=genus_sp), inherit.aes = F)+
  geom_point(size=10, shape=21, colour="black") +
  annotate("text", x = 2.25, y = 96, label = expression(paste(italic("Taygetis sp."), ' n=308')),size=9)+
  ylim(0,100)+ylab("Percetange abundance")+
  scale_fill_gradient(low = "#107a7c", high = "white", limits=c(0,100),
                      name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  scale_colour_gradient(low = "#107a7c", high = "white", limits=c(0,100), 
                        name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=27),
        strip.text.x = element_text(size = 18),
        axis.text =element_text(size = 20), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(colour = c("black","blue","black"), size = c(1,7,1)));taygetis.p 


italics.newline <- ~ atop(paste(italic("Chloreuptychia sp.")), paste("n=150"))

chloro.p <- ggplot(subset(sp.habitat.norm.plot.long , genus=="Chloreuptychia"&genus_sp%in%prim.sp), 
                   aes(x=habitat, y=norm.ab.perc, fill=abs(diff.to.prim), colour=abs(diff.to.prim)))+
  geom_line(aes(x=habitat, y=norm.ab.perc,group=genus_sp), inherit.aes = F)+
  geom_point(size=10, shape=21, colour="black") +
  #stat_summary(aes(y = norm.ab.perc,group=1), fun=mean, colour="red", geom="line",group=1)+
  annotate("text", x = 2.25, y = 92, label = italics.newline, size=9)+
  ylim(0,100)+ylab("Percetange abundance")+
  scale_fill_gradient(low = "#107a7c", high = "white", limits=c(0,100),
                      name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  scale_colour_gradient(low = "#107a7c", high = "white", limits=c(0,100), 
                        name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=27),
        strip.text.x = element_text(size = 18),
        axis.text =element_text(size = 20), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(colour = c("black","blue","black"), size = c(1,7,1)));chloro.p 




memphis.p <- ggplot(subset(sp.habitat.norm.plot.long , genus=="Pyrrhogyra"&genus_sp%in%prim.sp), 
                 aes(x=habitat, y=norm.ab.perc, fill=abs(diff.to.prim), colour=abs(diff.to.prim)))+
  geom_line(aes(x=habitat, y=norm.ab.perc,group=genus_sp), inherit.aes = F)+
  geom_point(size=10, shape=21, colour="black") +
  annotate("text", x = 2, y = 96, label = expression(paste(italic("Pyrrhogyra sp."), ' n=35')),size=9)+
  ylim(0,100)+ylab("Percetange abundance")+
  scale_fill_gradient(low = "#107a7c", high = "white", limits=c(0,100),
                      name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  scale_colour_gradient(low = "#107a7c", high = "white", limits=c(0,100), 
                        name="Difference between\n% of individuals\nfound in each habitat\nand % found in primary")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=27),
        strip.text.x = element_text(size = 18),
        axis.text =element_text(size = 20), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_text(size=22, face = c("plain", "bold","plain")),
        axis.ticks.x = element_line(colour = c("black","blue","black"), size = c(1,7,1)));memphis.p


plot_grid(taygetis.p, chloro.p, memphis.p, ncol=1, rel_heights = c(1,1,1.15))

# ggsave("plots/fig4/sp.plots.3.wide.png", height = 15, width = 7)

###### phylogenetic signal ######
# does preference for one or the other strategy show a phylogenetic signal??
# http://rfunctions.blogspot.fi/2014/02/measuring-phylogenetic-signal-in-r.html
# First, you need to define which trait you want to test and give names to each value according to species

################ primary species ################ 
#### mean change ####
sp.habitat.df.loss$mean.change <- (abs(sp.habitat.df.loss$Reserve)+abs(sp.habitat.df.loss$Logged))/2
trait <-data.frame(mean.change=sp.habitat.df.loss[,c("mean.change")]); rownames(trait)<- rownames(sp.habitat.df.loss); head(trait)
phylo <- as.phylo(tree.sophie.genus.sp.prim); plot(phylo)

# change tip names
rownames(trait)<-str_replace_all(str_replace_all(rownames(trait), "\\  ", "_"), "\\ ", "_")
df <- data.frame(genus.sp=phylo$tip.label, genus_sp=gsub(" ", "_",phylo$tip.label))
phylo <- rename_taxa(phylo,df, genus.sp, genus_sp )

phylo$tip.label; rownames(trait)

#rownames(trait)[11] <- "Archaeoprepona_demophoon"

#phylotraits <- NA
phylotraits <- phylo4d(phylo, trait); phylotraits

abouheif.test <- abouheif.moran(phylotraits,method="oriAbouheif"); abouheif.test
plot(abouheif.test)

#### riparian change ####
sp.habitat.df.loss$reserve.change <- abs(sp.habitat.df.loss$Reserve)
trait <-data.frame(reserve.change=sp.habitat.df.loss[,c("reserve.change")]); rownames(trait)<- rownames(sp.habitat.df.loss); head(trait)
phylo <- as.phylo(tree.genus.sp.prim); plot(phylo)
phylo <- as.phylo(tree.sophie.genus.sp.prim); plot(phylo)

# change tip names
rownames(trait)<-str_replace_all(str_replace_all(rownames(trait), "\\  ", "_"), "\\ ", "_")
df <- data.frame(genus.sp=phylo$tip.label, genus_sp=gsub(" ", "_",phylo$tip.label))
phylo <- rename_taxa(phylo,df, genus.sp, genus_sp )

phylo$tip.label; rownames(trait)
phylotraits <- phylo4d(phylo, trait); phylotraits

abouheif.test <- abouheif.moran(phylotraits,method="oriAbouheif"); abouheif.test
plot(abouheif.test)

#### logged only change ####
sp.habitat.df.loss$logged.change <- abs(sp.habitat.df.loss$Logged)
trait <-data.frame(logged.change=sp.habitat.df.loss[,c("logged.change")]); rownames(trait)<- rownames(sp.habitat.df.loss); head(trait)
phylo <- as.phylo(tree.sophie.genus.sp.prim); plot(phylo)

# change tip names
rownames(trait)<-str_replace_all(str_replace_all(rownames(trait), "\\  ", "_"), "\\ ", "_")
df <- data.frame(genus.sp=phylo$tip.label, genus_sp=gsub(" ", "_",phylo$tip.label))
phylo <- rename_taxa(phylo,df, genus.sp, genus_sp )

phylo$tip.label; rownames(trait)
phylotraits <- phylo4d(phylo, trait); phylotraits

abouheif.test <- abouheif.moran(phylotraits,method="oriAbouheif"); abouheif.test
plot(abouheif.test)

################ all species ################ 
tree.genus.sp.all 
tree.sophie.genus.sp.all

#### mean change ####
sp.habitat.df.loss$mean.change <- (abs(sp.habitat.df.loss$Reserve)+abs(sp.habitat.df.loss$Logged))/2
trait <-data.frame(mean.change=sp.habitat.df.loss[,c("mean.change")]); rownames(trait)<- rownames(sp.habitat.df.loss); head(trait)
phylo <- as.phylo(tree.sophie.genus.sp.all); plot(phylo)

# change tip names
rownames(trait)<-str_replace_all(str_replace_all(rownames(trait), "\\  ", "_"), "\\ ", "_")
df <- data.frame(genus.sp=phylo$tip.label, genus_sp=gsub(" ", "_",phylo$tip.label))
phylo <- rename_taxa(phylo,df, genus.sp, genus_sp )

phylo$tip.label; rownames(trait)
#rownames(trait)[11] <- "Archaeoprepona_demophoon"

#phylotraits <- NA
phylotraits <- phylo4d(phylo, trait); phylotraits

abouheif.test <- abouheif.moran(phylotraits,method="oriAbouheif"); abouheif.test
plot(abouheif.test)

#### riparian change ####
sp.habitat.df.loss$reserve.change <- abs(sp.habitat.df.loss$Reserve)
trait <-data.frame(reserve.change=sp.habitat.df.loss[,c("reserve.change")]); rownames(trait)<- rownames(sp.habitat.df.loss); head(trait)
phylo <- as.phylo(tree.sophie.genus.sp.all); plot(phylo)

# change tip names
rownames(trait)<-str_replace_all(str_replace_all(rownames(trait), "\\  ", "_"), "\\ ", "_")
df <- data.frame(genus.sp=phylo$tip.label, genus_sp=gsub(" ", "_",phylo$tip.label))
phylo <- rename_taxa(phylo,df, genus.sp, genus_sp )

phylo$tip.label; rownames(trait)
phylotraits <- phylo4d(phylo, trait); phylotraits

abouheif.test <- abouheif.moran(phylotraits,method="oriAbouheif"); abouheif.test
plot(abouheif.test)

#### logged change - significant ####
sp.habitat.df.loss$logged.change <- abs(sp.habitat.df.loss$Logged)
trait <-data.frame(logged.change=sp.habitat.df.loss[,c("logged.change")]); rownames(trait)<- rownames(sp.habitat.df.loss); head(trait)
phylo <- as.phylo(tree.sophie.genus.sp.all); plot(phylo)

# change tip names
rownames(trait)<-str_replace_all(str_replace_all(rownames(trait), "\\  ", "_"), "\\ ", "_")
df <- data.frame(genus.sp=phylo$tip.label, genus_sp=gsub(" ", "_",phylo$tip.label))
phylo <- rename_taxa(phylo,df, genus.sp, genus_sp )

phylo$tip.label; rownames(trait)
phylotraits <- phylo4d(phylo, trait); phylotraits

abouheif.test <- abouheif.moran(phylotraits,method="oriAbouheif"); abouheif.test
plot(abouheif.test)

################################################# S.I. #######################################################
########################## SOM fig S1 sp accum cruve ############ 
#### habitat level accumulation curves
#need to change df to aformat that inext can read. 
#species as rownames, groupsas columns titles and abundance data
names(sp.point)
sp.mat <- t(sp.habitat[, -c(1:3)]); sp.mat; colnames(sp.mat) <- t(sp.habitat[,2])

#### sp accum curve logged forest and prim comb #### 
sp.habitat$habitat<- c("PRI", "UPA.Int", "PRI","UPA.Str")
sp.mat.log <- t(subset(sp.habitat, substr(habitat, 0,3)=="UPA")[, -c(1:3)]); sp.mat.log; colnames(sp.mat.log) <- t(subset(sp.habitat, substr(habitat, 0,3)=="UPA")[,2])

sp.mat.pri.log <- t(subset(sp.habitat)[, -c(1:3)]); sp.mat.pri.log; colnames(sp.mat.pri.log) <- t(subset(sp.habitat)[,2]);sp.mat.pri.log

sp.mat.pri.log <- cbind(sp.mat.pri.log, sp.mat.pri.comb)
sp.mat.pri.log <-sp.mat.pri.log[,-c(1,3)];sp.mat.pri.log
colnames(sp.mat.pri.log)[3]<- "Prim"
head(sp.mat.pri.log)

#do the extrapolation with inext function
out.log<-iNEXT(sp.mat.pri.log, q=0, datatype="abundance", endpoint=1166*2)

# plot intrapolation (actual sp rich) and extrapolation with ggiNEXT
p.log <-ggiNEXT(out.log, type=1, color.var="site") +
  scale_color_manual(labels=c("Primary", "Interior logged","Riparian reserve"),
                     values=c( "#107a7c","#CC79A7", "#56B4E9")) +
  scale_fill_manual(labels=c("Primary", "Interior logged","Riparian reserve"),
                    values=c("#107a7c", "#CC79A7", "#56B4E9")) +
  scale_shape_manual(labels=c("Primary", "Interior logged","Riparian reserve"),values=c(15,17,19)) +
  theme(plot.background = element_blank(),legend.position="bottom" ,axis.text = element_text(size = 12),axis.title.x = element_text(),
        axis.title.y = element_text(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background=element_rect(colour="black", fill="white"), panel.spacing = element_blank()); p.log

# ggsave("plots/SI/fig.s1.png", width = 6, height = 6)

########################## SOM Fig. S2 forest structure by habitat type ############ 
#### habitat level accumulation curves

veg.point$habitat.location

veg.point$habitat.loc.3 <- if_else(substr(veg.point$site,0,3)=="PRI", "primary",veg.point$habitat.location ); veg.point$habitat.loc.3

veg.point$habitat.loc.3 <- factor(veg.point$habitat.loc.3, levels=c("primary", "UPA.Stream", "UPA.Interior"), labels=c("Primary", "Riparian\nreserve", "Interior\nlogged"))

names(veg.point)
my_comparisons1 <- list(c("Primary", "Interior\nlogged"), c("Primary", "Riparian\nreserve"), c("Interior\nlogged", "Riparian\nreserve") )

p1 <- ggplot(aes(x=habitat.loc.3, y=diam.tree.mean, fill=habitat.loc.3), data=veg.point)+
  geom_boxplot()+  
  stat_compare_means(method = "t.test", comparisons = my_comparisons1,label = "p.signif",
                                      tip.length = 0.0, vjust =0 , step.increase = 0.05, size=5) +
  scale_fill_manual(name="Forest type",  values=c("#107a7c" , "#56B4E9","#CC79A7"))+
  ylab("Mean tree diameter at breast height (cm)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 16),
        axis.text =element_text(size = 16), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_text(size=17, face = c("plain", "bold","plain")),
        axis.ticks.x = element_line(colour = c("black","blue","black"), size = c(1,7,1))); p1

p2 <- ggplot(aes(x=habitat.loc.3, y=mean.canopy, fill=habitat.loc.3), data=veg.point)+
  geom_boxplot()+  
  stat_compare_means(method = "t.test", comparisons = my_comparisons1,label = "p.signif",
                     tip.length = 0.0, vjust =0 , step.increase = 0.05, size=5) +
  scale_fill_manual(name="Forest type",  values=c("#107a7c" , "#56B4E9","#CC79A7"))+
  ylab("Mean canopy height (cm)")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 16),
        axis.text =element_text(size = 16), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_text(size=17, face = c("plain", "bold","plain")),
        axis.ticks.x = element_line(colour = c("black","blue","black"), size = c(1,7,1))); p2

p3 <- ggplot(aes(x=habitat.loc.3, y=canopy.cover, fill=habitat.loc.3), data=veg.point)+
  geom_boxplot()+  
  stat_compare_means(method = "t.test", comparisons = my_comparisons1,label = "p.signif",
                     tip.length = 0.0, vjust =0 , step.increase = 0.05, size=5) +
  scale_fill_manual(name="Forest type",  values=c("#107a7c" , "#56B4E9","#CC79A7"))+
  ylab("Canopy cover")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        strip.text.x = element_text(size = 16),
        axis.text =element_text(size = 16), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        #axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.text.x = element_text(size=17, face = c("plain", "bold","plain")),
        axis.ticks.x = element_line(colour = c("black","blue","black"), size = c(1,7,1))); p3


plot_grid(p1,p2, p3, labels = "AUTO", ncol=3)

# ggsave("plots/SI/forest.structure.png", width = 12, height = 6)

########################## SOM Note S1 recapture analyses ##########################
sp.recap <- read.csv("1.data/Butterflies.recaptured.csv", h = TRUE)
head(sp.recap)
head(sp)

library(stringi)
sp.recap$butt_code.original <- as.numeric(stri_extract_last(sp.recap$butt_code,regex="\\d+" ))

sp.recap$point_code<-gsub("PRI", "PRIM", sp.recap$point_code)
sp.recap$point_code.original <- sp$point_id[match(sp.recap$butt_code.original, sp$butt_code)]; sp.recap$point_code.original
sp.recap$site<-gsub("PRI", "PRIM", sp.recap$site)
sp.recap$site.original <- sp$site_id[match(sp.recap$site, sp$site_id)]; sp.recap$site.original

sp.recap$strip_id <- substr(sp.recap$point_code,0,7)
sp.recap$strip.original <- sp$strip_id[match(sp.recap$butt_code.original, sp$butt_code)]; sp.recap$strip.original


sp.recap$stream

length(unique(sp$site_id))
length(unique(sp$strip_id))

length(sp.recap$butt_code.original)
length(unique(sp.recap$butt_code.original))
uniq.recap <-unique(sp.recap$butt_code.original)
nrow(subset(sp.recap, butt_code.original %in% uniq.recap ))

recap.summ <- summarise(group_by(sp.recap,butt_code.original ),
                        n=n()); recap.summ
nrow(subset(recap.summ, n==1))
nrow(subset(recap.summ, n==2))
nrow(subset(recap.summ, n==3))
nrow(subset(recap.summ, n==4))

(nrow(subset(sp.recap, point_code==point_code.original))/ nrow(sp.recap) )*100 # 64% of recaptures were on same point 
(nrow(subset(sp.recap, strip_id==strip.original))/ nrow(sp.recap) )*100 # 100% of recaptures were within same site (never across transects)
(nrow(subset(sp.recap, site==site.original))/ nrow(sp.recap) )*100 # 100% of recaptures were within same site (never across transects)

sp.recap.diff.transect <- subset(sp.recap, strip_id!=strip.original)
length(unique(sp.recap.diff.transect$butt_code.original))

sp.recap.diff.transect$species <- sp$genus_sp[match(sp.recap.diff.transect$butt_code.original, sp$butt_code)]
sp.recap.diff.transect$genus <- sp$genus[match(sp.recap.diff.transect$butt_code.original, sp$butt_code)]

sp.recap.diff.transect$tribe <- sp$tribe[match(sp.recap.diff.transect$butt_code.original, sp$butt_code)]

summarise(group_by(unique(sp.recap.diff.transect[, c("butt_code.original", "genus", "species", "tribe")]), genus),
          n=n())

setwd("/Users/labadmin2/Documents/infant_gill/ProperLabels")
require(ggplot2)
require(ComplexHeatmap)
require(reshape2)
library(RColorBrewer)
require(DESeq2)
require(gridExtra)
require(Rmisc)
require(vegan)
require(circlize)
require(hclust)
require(NbClust)
require(coda.base)
require(usedist)

#Create color palette to be used for plotting
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Load the count matrixes for phylum, genus, and species
phyla <- readRDS("RefSeq2_p.RDS")
genus <- readRDS("RefSeq2_g.RDS")
species <- readRDS("RefSeq2_s.RDS")

#delta/epsilon subdivisions is a part of Proteobacteria, so we're going with that.
phyla["Proteobacteria",] <- phyla["Proteobacteria",] + phyla["delta/epsilon subdivisions",]
phyla <- phyla[rownames(phyla) != "delta/epsilon subdivisions",]
#Require more than 1000 reads in a sample to be analyzed
phyla <- phyla[,apply(phyla,2,sum)>10000]
phyla <- phyla[apply(phyla,1,sum) > 10,]
genus <- genus[,apply(genus,2,sum)>10000]
genus <- genus[apply(genus,1,sum) > 10,]
species <- species[,apply(species,2,sum)>10000]
species <- species[apply(species,1,sum) > 10,]
colnames(genus) <- gsub("S|S0","X",colnames(genus))
colnames(phyla) <- gsub("S|S0","X",colnames(phyla))
colnames(species) <- gsub("S|S0","X",colnames(species))
colnames(genus)[grep("DNA",colnames(genus))] <- c("DNA1","DNA2","DNA3","DNA4","DNA5")
colnames(phyla)[grep("DNA",colnames(phyla))] <- c("DNA1","DNA2","DNA3","DNA4","DNA5")
colnames(species)[grep("DNA",colnames(species))] <- c("DNA1","DNA2","DNA3","DNA4","DNA5")
#Create a counts-per-million matrix for each taxonomic level
CPMphyla <- t(t(phyla)/apply(phyla,2,sum)*1000000)
CPMgenus <- t(t(genus)/apply(genus,2,sum)*1000000)
CPMspecies <- t(t(species)/apply(species,2,sum)*1000000)
map <- readRDS("../InfantMappingData2.RDS")
map$PrePostSympt[is.na(map$PrePostSympt)] <- "none"
map$PrePostSympt[map$PrePostSympt=="none"] <- "Healthy"
map$PrePostSympt[map$Symptoms=="Yes" & !is.na(map$Symptoms)] <- "Active"

fullGenus <- readRDS("ggGenusWithClusters.RDS")

BetaDivMap <- dist(t(CPMgenus[apply(CPMgenus[,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000]],1,sum)>0,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000]])+0.1,method="aitchison")
BetaDivMap <- dist_setNames(BetaDivMap, colnames(CPMgenus[apply(CPMgenus[,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000]],1,sum)>0,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000]]))
healthyBetaDivMap <- dist(t(CPMgenus[apply(CPMgenus[,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000 & map$InfectionStatus=="Healthy"]],1,sum)>0,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000 & map$InfectionStatus=="Healthy"]])+0.1,method="aitchison")
healthyBetaDivMap <- dist_setNames(healthyBetaDivMap,colnames(CPMgenus[apply(CPMgenus[,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000 & map$InfectionStatus=="Healthy"]],1,sum)>0,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000 & map$InfectionStatus=="Healthy"]]))
sickBetaDivMap <- dist(t(CPMgenus[apply(CPMgenus[,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000 & map$InfectionStatus=="RespiratoryIllness"]],1,sum)>0,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000 & map$InfectionStatus=="RespiratoryIllness"]])+0.1,method="aitchison")
sickBetaDivMap <- dist_setNames(sickBetaDivMap,colnames(CPMgenus[apply(CPMgenus[,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000 & map$InfectionStatus=="RespiratoryIllness"]],1,sum)>0,map$Sample[map$MothChild=="Infant" & map$Sample %in% colnames(CPMgenus) & map$Age < 1000 & map$InfectionStatus=="RespiratoryIllness"]]))

MyDendroHealthy <- hclust(healthyBetaDivMap, method="ward.D")
#Silhouette index maxes out at 5 clusters
plot(3:21,NbClust(diss=healthyBetaDivMap, distance=NULL, method="ward.D", min.nc=3, max.nc=21, index="silhouette")$All.index)
#Frey index maxes out at 15 clusters (huge outlier)
plot(3:21,NbClust(diss=healthyBetaDivMap, distance=NULL, method="ward.D", min.nc=3, max.nc=21, index="frey")$All.index)

#Should label 5/15 clusters for healthy infants only
healthyClusters5 <- cutree(MyDendroHealthy, k=5)
healthyClusters15 <- cutree(MyDendroHealthy, k=15)
col_fun = colorRamp2(c(0,60,120),c("white","yellow","black"))

healthyTopAnno <- HeatmapAnnotation(Age=map[colnames(as.matrix(healthyBetaDivMap)),"Age"],
                                    #HIV_Exposure=map$HIVStatus[map$Sample %in% colnames(as.matrix(healthyBetaDivMap))],
                                    Cluster=healthyClusters5,
                                    MoreCluster=healthyClusters15,
                                    col=list(Age=col_fun,
                                             #HIV_Exposure=c("Control"="blue","HIV"="red"),
                                             Cluster=c("1"=col_vector[1],"2"=col_vector[2],"3"=col_vector[3],"4"=col_vector[4],"5"=col_vector[5]),
                                             MoreCluster=c("1"=col_vector[1],"2"=col_vector[2],"3"=col_vector[3],"4"=col_vector[4],"5"=col_vector[5],
                                                           "6"=col_vector[6],"7"=col_vector[7],"8"=col_vector[8],"9"=col_vector[9],"10"=col_vector[10],
                                                           "11"=col_vector[11],"12"=col_vector[12],"13"=col_vector[13],"14"=col_vector[14],"15"=col_vector[15])))
Heatmap(as.matrix(healthyBetaDivMap), top_annotation = healthyTopAnno, show_row_names = FALSE,
        show_column_names = FALSE, name = "Aitchison distance", cluster_rows = MyDendroHealthy, cluster_columns = MyDendroHealthy,
        column_title = "Aitchison distance matrix", column_title_gp = gpar(fontsize = 20, fontface = "bold"))

fullGenus$HC5 <- 0
for(i in 1:5){
  fullGenus$HC5[fullGenus$Sample %in% names(healthyClusters5)[healthyClusters5==i]] <- i
}
fullGenus$HC15 <- 0
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 1]] <- "1a"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 5]] <- "1b"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 7]] <- "1c"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 14]] <- "1d"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 2]] <- "2"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 12]] <- "3b"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 13]] <- "3c"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 15]] <- "3d"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 3]] <- "3a"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 4]] <- "4a"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 8]] <- "4b"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 9]] <- "4c"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 10]] <- "5a"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 11]] <- "5b"
fullGenus$HC15[fullGenus$Sample %in% names(healthyClusters15)[healthyClusters15 == 6]] <- "5c"


require(ape)
healthyPCoA <- pcoa(healthyBetaDivMap)
plotHealthyPCoA <- data.frame(healthyPCoA$vectors[,1:6])
colnames(plotHealthyPCoA) = c("PC_1","PC_2","PC_3","PC_4","PC_5","PC_6")
plotHealthyPCoA$Subject <- factor(map[rownames(healthyPCoA$vectors),"Subject"])
plotHealthyPCoA$Age <- map[rownames(healthyPCoA$vectors),"Age"]
# head(healthyPCoA$values)
ggplot(plotHealthyPCoA, aes(x=PC_1, y=PC_2, fill=Age)) +
  geom_point(size=5, shape=21, color="black") +
  theme_linedraw() +
  scale_fill_gradient2(name="Age",low="white",mid="yellow",high="black",midpoint=60) +
  theme(axis.title = element_text(face="bold", colour="#000000", size=18)) +
  theme(legend.text = element_text(colour="grey35", size = 18, face = "bold.italic")) +
  theme(legend.title = element_text(colour="grey20", size = 22, face = "bold")) + 
  xlab("PC1 (RCE:0.06)") +
  ylab("PC2 (RCE:0.05)") +
  theme(axis.text.x  = element_text(size=20, angle = 90, vjust = 0.5)) +
  theme(axis.text.y  = element_text(size=20)) +
  ggtitle("PCoA (on  Bray-Curtis distances)\nof healthy infants") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))



MyDendro <- hclust(BetaDivMap, method="ward.D")
#plot(3:21,NbClust(diss=BetaDivMap, distance=NULL, method="ward.D", min.nc=3, max.nc=21, index="silhouette")$All.index)
#plot(3:21,NbClust(diss=BetaDivMap, distance=NULL, method="ward.D", min.nc=3, max.nc=21, index="frey")$All.index)
# NbClust(diss=BetaDivMap, distance=NULL, method="average", min.nc=5, max.nc=15, index="silhouette")
# NbClust(diss=BetaDivMap, distance=NULL, method="average", min.nc=2, max.nc=15, index="frey")
basicClusters <- cutree(MyDendro, k=6)
basicClusters2 <- cutree(MyDendro, k=11)
TopAnno <- HeatmapAnnotation(Symptoms=map$PrePostSympt[map$Sample %in% colnames(as.matrix(BetaDivMap))],
                             Age=map$Age[map$Sample %in% colnames(as.matrix(BetaDivMap))],
                             #HIV_Exposure=map$HIVStatus[map$Sample %in% colnames(as.matrix(BetaDivMap))],
                             Cluster=basicClusters,
                             MoreCluster=basicClusters2,
                             col=list(Symptoms=c("Pre"="goldenrod","Post"="darkgreen","Healthy"="lightblue", "Active"="darkred"),
                                      Age=col_fun,
                                      #HIV_Exposure=c("Control"="blue","HIV"="red"),
                                      Cluster=c("1"=col_vector[1],"2"=col_vector[2],"3"=col_vector[3],"4"=col_vector[4],"5"=col_vector[5],"6"=col_vector[6]),
                                      MoreCluster=c("1"=col_vector[1],"2"=col_vector[2],"3"=col_vector[3],"4"=col_vector[4],"5"=col_vector[5],"6"=col_vector[6],"7"=col_vector[7],"8"=col_vector[8],"9"=col_vector[9],"10"=col_vector[10],
                                                    "11"=col_vector[11])))
Heatmap(as.matrix(BetaDivMap), top_annotation = TopAnno, show_row_names = FALSE,
        show_column_names = FALSE, name = "Aitchison", cluster_rows = MyDendro, cluster_columns = MyDendro)
```

#### PCoA {.tabset}

##### Age

```{r, echo=FALSE, warning=FALSE, message=FALSE}
require(ape)
myPCoA <- pcoa(BetaDivMap)
plotAllPCoA <- data.frame(myPCoA$vectors[,1:6])
colnames(plotAllPCoA) = c("PC_1","PC_2","PC_3","PC_4","PC_5","PC_6")
plotAllPCoA$Subject <- factor(map[rownames(myPCoA$vectors),"Subject"])
plotAllPCoA$Age <- map[rownames(myPCoA$vectors),"Age"]
plotAllPCoA$LRTI <- map[rownames(myPCoA$vectors),"PrePostSympt"]
# head(myPCoA$values)
ggplot(plotAllPCoA, aes(x=PC_1, y=PC_2, fill=Age)) +
  geom_point(size=5, shape=21, color="black") +
  theme_linedraw() +
  scale_fill_gradient2(name="Age",low="white",mid="yellow",high="black",midpoint=60) +
  theme(axis.title = element_text(face="bold", colour="#000000", size=18)) +
  theme(legend.text = element_text(colour="grey35", size = 18, face = "bold.italic")) +
  theme(legend.title = element_text(colour="grey20", size = 22, face = "bold")) + 
  xlab("PC1 (RCE:0.05)") +
  ylab("PC2 (RCE:0.04)") +
  theme(axis.text.x  = element_text(size=20, angle = 90, vjust = 0.5)) +
  theme(axis.text.y  = element_text(size=20)) +
  ggtitle("PCoA (on  Bray-Curtis distances)\nof all infant samples") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))
ggplot(plotAllPCoA, aes(x=PC_3, y=PC_4, fill=Age)) +
  geom_point(size=5, shape=21, color="black") +
  theme_linedraw() +
  scale_fill_gradient2(name="Age",low="white",mid="yellow",high="black",midpoint=60) +
  theme(axis.title = element_text(face="bold", colour="#000000", size=18)) +
  theme(legend.text = element_text(colour="grey35", size = 18, face = "bold.italic")) +
  theme(legend.title = element_text(colour="grey20", size = 22, face = "bold")) + 
  xlab("PC3 (RCE:0.03)") +
  ylab("PC4 (RCE:0.02)") +
  theme(axis.text.x  = element_text(size=20, angle = 90, vjust = 0.5)) +
  theme(axis.text.y  = element_text(size=20)) +
  ggtitle("PCoA (on  Bray-Curtis distances)\nof all infant samples") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))
ggplot(plotAllPCoA, aes(x=PC_5, y=PC_6, fill=Age)) +
  geom_point(size=5, shape=21, color="black") +
  theme_linedraw() +
  scale_fill_gradient2(name="Age",low="white",mid="yellow",high="black",midpoint=60) +
  theme(axis.title = element_text(face="bold", colour="#000000", size=18)) +
  theme(legend.text = element_text(colour="grey35", size = 18, face = "bold.italic")) +
  theme(legend.title = element_text(colour="grey20", size = 22, face = "bold")) + 
  xlab("PC5 (RCE:0.02)") +
  ylab("PC6 (RCE:0.01)") +
  theme(axis.text.x  = element_text(size=20, angle = 90, vjust = 0.5)) +
  theme(axis.text.y  = element_text(size=20)) +
  ggtitle("PCoA (on  Bray-Curtis distances)\nof all infant samples") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))
```

##### LRTI status

```{r, echo=FALSE, warning=FALSE, message=FALSE}
ggplot(plotAllPCoA, aes(x=PC_1, y=PC_2, fill=factor(LRTI, levels=c("Healthy","Pre","Active","Post")))) +
  geom_point(size=5, shape=21, color="black") +
  theme_linedraw() +
  scale_fill_manual(values=c("lightblue","goldenrod","darkred","darkgreen")) +
  guides(fill = guide_legend(ncol = 1, bycol=TRUE, title="LRTI status")) + 
  theme(axis.title = element_text(face="bold", colour="#000000", size=18)) +
  theme(legend.text = element_text(colour="grey35", size = 18, face = "bold.italic")) +
  theme(legend.title = element_text(colour="grey20", size = 22, face = "bold")) + 
  xlab("PC1 (RCE:0.05)") +
  ylab("PC2 (RCE:0.04)") +
  theme(axis.text.x  = element_text(size=20, angle = 90, vjust = 0.5)) +
  theme(axis.text.y  = element_text(size=20)) +
  ggtitle("PCoA (on  Bray-Curtis distances)\nof all infant samples") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))
ggplot(plotAllPCoA, aes(x=PC_3, y=PC_4, fill=factor(LRTI, levels=c("Healthy","Pre","Active","Post")))) +
  geom_point(size=5, shape=21, color="black") +
  theme_linedraw() +
  scale_fill_manual(values=c("lightblue","goldenrod","darkred","darkgreen")) +
  guides(fill = guide_legend(ncol = 1, bycol=TRUE, title="LRTI status")) + 
  theme(axis.title = element_text(face="bold", colour="#000000", size=18)) +
  theme(legend.text = element_text(colour="grey35", size = 18, face = "bold.italic")) +
  theme(legend.title = element_text(colour="grey20", size = 22, face = "bold")) + 
  xlab("PC3 (RCE:0.03)") +
  ylab("PC4 (RCE:0.02)") +
  theme(axis.text.x  = element_text(size=20, angle = 90, vjust = 0.5)) +
  theme(axis.text.y  = element_text(size=20)) +
  ggtitle("PCoA (on  Bray-Curtis distances)\nof all infant samples") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))
ggplot(plotAllPCoA, aes(x=PC_5, y=PC_6, fill=factor(LRTI, levels=c("Healthy","Pre","Active","Post")))) +
  geom_point(size=5, shape=21, color="black") +
  theme_linedraw() +
  scale_fill_manual(values=c("lightblue","goldenrod","darkred","darkgreen")) +
  guides(fill = guide_legend(ncol = 1, bycol=TRUE, title="LRTI status")) + 
  theme(axis.title = element_text(face="bold", colour="#000000", size=18)) +
  theme(legend.text = element_text(colour="grey35", size = 18, face = "bold.italic")) +
  theme(legend.title = element_text(colour="grey20", size = 22, face = "bold")) + 
  xlab("PC5 (RCE:0.02)") +
  ylab("PC6 (RCE:0.01)") +
  theme(axis.text.x  = element_text(size=20, angle = 90, vjust = 0.5)) +
  theme(axis.text.y  = element_text(size=20)) +
  ggtitle("PCoA (on  Bray-Curtis distances)\nof all infant samples") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))


miniCPMgg <- fullGenus[fullGenus$Genus=="Other",]

miniCPMgg$logValue <- log(miniCPMgg$Value*10000+1)


ggplot(miniCPMgg[miniCPMgg$Age < 1000,],
       aes(x=factor(PrePost,levels=c("Healthy","Pre","Active","Post")), y=logValue,
           fill=factor(PrePost,levels=c("Healthy","Pre","Active","Post")))) +
  stat_boxplot(geom ='errorbar',width = 0.3, size=1) +
  theme(axis.line=element_line(color="black",size=0.5)) +
  geom_boxplot(lwd=1, color="#000000") +
#  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_fill_manual(values=c("lightblue","goldenrod","darkred","darkgreen")) +
  guides(fill = guide_legend(ncol = 1,bycol=TRUE,title="LRTI status")) + 
  theme(axis.title = element_text(face="bold", colour="#000000", size=18)) +
  xlab("LRTI status") +
  ylab("log(CPM)") +
  theme(axis.text.x  = element_text(size=14, face="bold", angle = 0, hjust = 0.5)) +
  theme(axis.text.y  = element_text(size=14)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,16)) +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill="gray95")) +
  ggtitle("Relative abundance of rare genera") +
  theme(plot.title = element_text(hjust = 0.5,face="bold",size=20))

t.test(miniCPMgg$logValue[miniCPMgg$Age<1000 & miniCPMgg$PrePost=="RespiratoryIllness"],
       miniCPMgg$logValue[miniCPMgg$Age<1000 & miniCPMgg$InfectionStatus=="Healthy"])


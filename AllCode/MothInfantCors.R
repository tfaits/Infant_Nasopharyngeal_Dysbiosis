setwd("/Users/labadmin2/Documents/infant_gill/ProperLabels")
require(ggplot2)
require(ComplexHeatmap)
require(reshape2)
library(RColorBrewer)
require(DESeq2)
require(gridExtra)
require(Rmisc)
require(vegan)
require(fossil)
require(circlize)
require(knitr)
require(lme4)
require(car)
phyla <- readRDS("phyla.RDS")
genus <- readRDS("genus.RDS")
species <- readRDS("species.RDS")
CPMggPhyla <- readRDS("CPMggPhyla.RDS")
CPMggGenus <- readRDS("CPMggGenus.RDS")
CPMggSpecies <- readRDS("CPMggSpecies.RDS")
CPMphyla <- readRDS("CPMphyla.RDS")
CPMgenus <- readRDS("CPMgenus.RDS")
CPMspecies <- readRDS("CPMspecies.RDS")
map <- readRDS("map.RDS")
#Create color palette to be used for plotting
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
sickInfants <- unique(map$Subject[map$MothChild=="Infant" & map$InfectionStatus=="RespiratoryIllness"])

FirstMoms <- c()
for(i in unique(map[map$MothChild=="Mother","Subject"])){
  tmp <- min(as.numeric(gsub("X","",map$Sample[map$Subject==i])))
  FirstMoms <- c(FirstMoms, paste("X",tmp,sep=""))
}

MothKids <- list()
for(Mother in FirstMoms){
  tmp <- map$Sample[map$Subject==gsub("-0","-1",map$Subject[map$Sample==Mother])]
  MothKids[[Mother]] <- tmp[order(as.numeric(gsub("X","",tmp)))]
}

#topGen <- c("Moraxella","Streptococcus","Haemophilus","Dolosigranulum","Bacillus","Paracoccus",
#            "Acinetobacter","Anaerobacillus","Pseudomonas","Staphylococcus","Corynebacterium",
#            "Delftia","Novosphingobium")

topGen <- c("Staphylococcus","Streptococcus","Corynebacterium","Acinetobacter",
            "Paracoccus","Delftia","Bacillus","Pseudomonas","Dolosigranulum","Moraxella")

MaternalInfantCor <- list()
for(genus in topGen){
  MaternalInfantCor[[genus]] <- data.frame(Mom=log(CPMgenus[genus,names(MothKids)]+1),
                                           Child=log(CPMgenus[genus,unlist(lapply(MothKids,function(x){x[1]}))]+1))
}

tmpCor <- data.frame(genus=c("None"),Spearman=c(0),pVal=c(1),stringsAsFactors = FALSE)
for(genus in topGen){
  tmp <- cor(MaternalInfantCor[[genus]][,1],MaternalInfantCor[[genus]][,2])
  randCor <- c()
  for(i in 1:10000){
    randCor <- c(randCor,cor(MaternalInfantCor[[genus]][sample(1:40,40,replace=TRUE),1],MaternalInfantCor[[genus]][sample(1:40,40,replace=TRUE),2]))
  }
  tmpCor <- rbind(tmpCor,c(genus,tmp,(sum(abs(randCor)>=abs(tmp))+1)/10000))
}
tmpCor <- tmpCor[-1,]
tmpCor$qVal <- p.adjust(tmpCor$pVal)

#P-value that the overall correlations are not 0:
t.test(as.numeric(tmpCor$Spearman))

## Do the same correlations, but looking at the mean value in the infants over time
MatInfCorMean <- list()
for(genus in topGen){
  MatInfCorMean[[genus]] <- data.frame(Mom=log(CPMgenus[genus,names(MothKids)]+1),
                                       Child=unlist(lapply(MothKids,function(x){log(mean(CPMgenus[genus,x])+1)})))
}
MeanCor <- data.frame(genus=c("None"),Spearman=c(0),pVal=c(1),stringsAsFactors = FALSE)
for(genus in topGen){
  tmp <- cor(MatInfCorMean[[genus]][,1],MatInfCorMean[[genus]][,2])
  randCor <- c()
  for(i in 1:10000){
    randCor <- c(randCor,cor(MatInfCorMean[[genus]][sample(1:40,40,replace=TRUE),1],MatInfCorMean[[genus]][sample(1:40,40,replace=TRUE),2]))
  }
  MeanCor <- rbind(MeanCor,c(genus,tmp,(sum(abs(randCor)>=abs(tmp))+1)/10000))
}
MeanCor <- MeanCor[-1,]
MeanCor$qVal <- p.adjust(MeanCor$pVal)


## Do the same correlations, but looking at the maximum value in the infants over time
MatInfCorMax <- list()
for(genus in topGen){
  MatInfCorMax[[genus]] <- data.frame(Mom=log(CPMgenus[genus,names(MothKids)]+1),
                                       Child=unlist(lapply(MothKids,function(x){log(max(CPMgenus[genus,x])+1)})))
}
MaxCor <- data.frame(genus=c("None"),Spearman=c(0),pVal=c(1),stringsAsFactors = FALSE)
for(genus in topGen){
  tmp <- cor(MatInfCorMax[[genus]][,1],MatInfCorMax[[genus]][,2])
  randCor <- c()
  for(i in 1:10000){
    randCor <- c(randCor,cor(MatInfCorMax[[genus]][sample(1:40,40,replace=TRUE),1],MatInfCorMax[[genus]][sample(1:40,40,replace=TRUE),2]))
  }
  MaxCor <- rbind(MaxCor,c(genus,tmp,(sum(abs(randCor)>=abs(tmp))+1)/10000))
}
MaxCor <- MaxCor[-1,]
MaxCor$qVal <- p.adjust(MaxCor$pVal)

BetaDivMap <- vegdist(t(CPMgenus[,map$Age < 1000 | map$Sample %in% FirstMoms]), method="bray")
MyDendro <- hclust(BetaDivMap, method="average")
numClusts <- 11
basicClusters <- cutree(MyDendro, k=numClusts)
#Creating some colors - not important and very much your own aesthetic choices
col_fun = colorRamp2(c(0,60,120),c("white","yellow","black"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ClusterColFun <- col_vector[1:numClusts]
TopAnno <- HeatmapAnnotation(HIVStatus=map[colnames(as.matrix(BetaDivMap)),"HIVStatus"],
                             Symptoms=map[colnames(as.matrix(BetaDivMap)),"PrePostSympt"],
                             Age=unlist(lapply(colnames(as.matrix(BetaDivMap)),function(x){if(map[x,"MothChild"]=="Infant"){map[x,"Age"]}else{NA}})),
                             Cluster=basicClusters,
                             MothChild=map[colnames(as.matrix(BetaDivMap)),"MothChild"],
                             col=list(HIVStatus=c("HIV"="red","Control"="blue"),
                                      Symptoms=c("Pre"="goldenrod","Post"="darkgreen",
                                                 "Active"="darkred",Healthy="lightblue"),
                                      Age=col_fun,
                                      Cluster=c("1"=col_vector[1],"2"=col_vector[2],"3"=col_vector[3],
                                                "4"=col_vector[4],"5"=col_vector[5],"6"=col_vector[6],
                                                "7"=col_vector[7],"8"=col_vector[8],"9"=col_vector[9],
                                                "10"=col_vector[10],"11"=col_vector[11]
                                      ),
                                      MothChild=c("Infant"="pink","Mother"="green")
                             )
)

Heatmap(as.matrix(BetaDivMap), top_annotation = TopAnno, show_row_names = FALSE,
        show_column_names = FALSE, name = "Bray-Curtis", cluster_rows = MyDendro, cluster_columns = MyDendro)

betaDivMap <- as.matrix(BetaDivMap)
#Boxplot of the BC distance from each mother to her own child vs each mother to all other children
motherToOwn <- unlist(lapply(names(MothKids),function(x){as.matrix(BetaDivMap)[x,MothKids[[x]][1]]}))
boxplot(motherToOwn,as.matrix(BetaDivMap)[names(MothKids),unlist(lapply(MothKids,function(x){x[1]}))])

#Boxplot for each mother - her distance to other kids vs her own:
#Also keep track of rank
similarityRank <- c()
shortDist <- c()
AllDist <- c()
for(mom in names(MothKids)){
  #myPlot <- boxplot(betaDivMap[mom,MothKids[[mom]][1]],betaDivMap[mom,unlist(lapply(names(MothKids)[names(MothKids)!=mom],function(x){MothKids[[x]][1]}))])
  #print(myPlot)
  similarityRank <- c(similarityRank, sum(betaDivMap[mom,unlist(lapply(names(MothKids)[names(MothKids)!=mom],function(x){MothKids[[x]][1]}))]<betaDivMap[mom,MothKids[[mom]][1]]))
  shortDist <- c(shortDist,betaDivMap[mom,MothKids[[mom]][1]])
  AllDist <- c(AllDist, betaDivMap[mom,unlist(lapply(names(MothKids)[names(MothKids)!=mom],function(x){MothKids[[x]][1]}))])
}

ks.test(x=similarityRank,y=sample(0:39,10000000,replace=TRUE))

ks.test(x=shortDist, y=AllDist)


###Differences between mothers
require(DESeq2)
mothsPhy <- DESeqDataSetFromMatrix(countData = round(phyla[,FirstMoms]) + 1,
                                                colData = map[FirstMoms,],
                                                                design = ~ HIVStatus + InfectionStatus)
mothsPhy <- DESeq(mothsPhy)
saveRDS(mothsPhy,"mothsPhy.RDS")
#mothsPhy <- readRDS("mothsPhy.RDS")
mothsResphy <- results(mothsPhy)
mothsResultphy <- as.matrix(mothsResphy[, c("baseMean","log2FoldChange","pvalue","padj")])
mothsResultphy[is.na(mothsResultphy[,"padj"]),"padj"] <- 1
mothsDifferencesphy <- mothsResultphy[mothsResultphy[,"baseMean"] > 5,]

mothsGen <- DESeqDataSetFromMatrix(countData = round(genus[apply(genus[,FirstMoms],1,function(x){sum(x>10)})>3,FirstMoms]) + 1,
                                   colData = map[FirstMoms,],
                                   design = ~ HIVStatus + InfectionStatus)
mothsGen <- DESeq(mothsGen)
#saveRDS(mothsGen,"mothsGen.RDS")
#mothsGen <- readRDS("mothsGen.RDS")
mothsResgen <- results(mothsGen)
mothsResultgen <- as.matrix(mothsResgen[, c("baseMean","log2FoldChange","pvalue","padj")])
mothsResultgen[is.na(mothsResultgen[,"padj"]),"padj"] <- 1
mothsDifferencesgen <- mothsResultgen[mothsResultgen[,"padj"] < 0.1,]
mothsDifferencesgen

plotGens <- function(taxon){
  boxplot((CPMgenus[taxon,FirstMoms[FirstMoms %in% map$Sample[map$InfectionStatus=="Healthy"]]]+1),
          (CPMgenus[taxon,FirstMoms[FirstMoms %in% map$Sample[map$InfectionStatus=="RespiratoryIllness"]]]+1),
          names=c("Healthy","SRTI"),
          main=paste(taxon,"abundance"))
}
plotGens2 <- function(taxon){
  boxplot((CPMgenus[taxon,FirstMoms[FirstMoms %in% map$Sample[map$InfectionStatus=="Healthy" & map$HIVStatus=="Control"]]]+1),
          (CPMgenus[taxon,FirstMoms[FirstMoms %in% map$Sample[map$InfectionStatus=="RespiratoryIllness" & map$HIVStatus=="Control"]]]+1),
          names=c("Healthy","SRTI"),
          main=paste(taxon,"abundance"))
}

mothsResultgen[c("Anaerobacillus","Bacillus","Blastococcus","Brachybacterium",
                 "Delftia","Dolosigranulum","Novosphingobium","Ochrobactrum",
                 "Ornithinimicrobium","Sphingomonas")]


mothsPhy <- DESeqDataSetFromMatrix(countData = round(phyla[,FirstMoms]) + 1,
                                   colData = map[FirstMoms,],
                                   design = ~ HIVStatus + InfectionStatus)
mothsPhy <- DESeq(mothsPhy)
saveRDS(mothsPhy,"mothsPhy.RDS")
#mothsPhy <- readRDS("mothsPhy.RDS")
mothsResphy <- results(mothsPhy)
mothsResultphy <- as.matrix(mothsResphy[, c("baseMean","log2FoldChange","pvalue","padj")])
mothsResultphy[is.na(mothsResultphy[,"padj"]),"padj"] <- 1
mothsDifferencesphy <- mothsResultphy[mothsResultphy[,"baseMean"] > 5,]



#Beta Diversity ordering?
binaryGenus <- CPMgenus
for(taxon in rownames(binaryGenus)){
  binaryGenus[taxon,] <- unlist(lapply(CPMgenus[taxon,],function(x){if(x>10){return(1)}else{return(0)}}))
}
jacDivMap <- vegdist(t(binaryGenus[,map$Age < 1000 | map$Sample %in% FirstMoms]), method="jaccard")
jacdivmap <- as.matrix(jacDivMap)

firstKids <- map$Sample[map$timepoint==0 & !is.na(map$timepoint)]
MomToKid <- unlist(lapply(names(MothKids),function(x){jacdivmap[x,MothKids[[x]][1]]}))
MomToOthers <- lapply(names(MothKids),function(x){jacdivmap[x,firstKids[firstKids!=MothKids[[x]][1]]]})

momDistRanks <- unlist(lapply(1:40,function(x){sum(MomToOthers[[x]] < MomToKid[x])}))

boxplot(MomToKid,unlist(MomToOthers))

#WorstCorrelation
names(MothKids)[momDistRanks>30]

#BestCorrelation
names(MothKids)[momDistRanks<10]

ks.test(momDistRanks,"punif",0,39)

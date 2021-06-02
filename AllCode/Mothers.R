

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
require(limma)
phyla <- readRDS("phyla.RDS")
genus <- readRDS("genus.RDS")
species <- readRDS("species.RDS")
CPMphyla <- t(t(phyla)/apply(phyla,2,sum))*1000000
CPMgenus <- t(t(genus)/apply(genus,2,sum))*1000000
CPMspecies <- t(t(species)/apply(species,2,sum))*1000000
CPMggPhyla <- readRDS("CPMggPhyla.RDS")
CPMggGenus <- readRDS("CPMggGenus.RDS")
CPMggSpecies <- readRDS("CPMggSpecies.RDS")
map <- readRDS("map.RDS")

firstMotherSamples <- c()
for(mother in unique(map$Subject[map$MothChild=="Mother"])){
  firstMotherSamples <- c(firstMotherSamples,map$Sample[map$Subject==mother][order(as.numeric(gsub("X","",map$Sample[map$Subject==mother])))[1]])
}

mothGen <-round(genus[,firstMotherSamples])
mothGen <- mothGen[apply(mothGen,1,sum) >= 100 & apply(mothGen,1,function(x){sum(x>0)}) > 3,]

motherGen <- DESeqDataSetFromMatrix(countData = round(mothGen),
                                  colData = map[colnames(mothGen),],
                                  design = ~ HIV + InfectionStatus)
#motherDE <- DESeq(motherGen)
#saveRDS(motherDE,"motherDE.RDS")
motherDE <- readRDS("motherDE.RDS")
motherRes <- results(motherDE)
motherResult <- as.matrix(motherRes[, c("baseMean","log2FoldChange","pvalue","padj")])
motherResult[is.na(motherResult[,"padj"]),"padj"] <- 1
motherResult[motherResult[,"baseMean"] > 500,]
motherDifferences <- motherResult[motherResult[,"padj"] < 0.1 & (motherResult[,"baseMean"] > 1000 | motherResult[,"baseMean"]*2^motherResult[,"log2FoldChange"] > 1000),]
kable(motherDifferences)

makeBoxPlot <- function(taxon){
  boxplot(CPMgenus[taxon,map$InfectionStatus=="Healthy" & map$Sample %in% firstMotherSamples]/10000,
          CPMgenus[taxon,map$InfectionStatus=="RespiratoryIllness" & map$Sample %in% firstMotherSamples]/10000,
          main=paste(taxon,"abundance in mothers"),names=c("Healthy","LRTI"))
}



### Try it with all mother samples
mothGen <-round(genus[,map$Sample[map$MothChild=="Mother"]])
mothGen <- mothGen[apply(mothGen,1,sum) >= 100 & apply(mothGen,1,function(x){sum(x>0)}) > 3,]

motherGen <- DESeqDataSetFromMatrix(countData = round(mothGen),
                                    colData = map[colnames(mothGen),],
                                    design = ~ HIV + InfectionStatus)
motherDE <- DESeq(motherGen)
motherRes <- results(motherDE)
motherResult <- as.matrix(motherRes[, c("baseMean","log2FoldChange","pvalue","padj")])
motherResult[is.na(motherResult[,"padj"]),"padj"] <- 1
motherResult[motherResult[,"baseMean"] > 500,]
motherDifferences <- motherResult[motherResult[,"padj"] < 0.1 & (motherResult[,"baseMean"] > 1000 | motherResult[,"baseMean"]*2^motherResult[,"log2FoldChange"] > 1000),]
kable(motherDifferences)

makeBoxPlot <- function(taxon){
  boxplot(CPMgenus[taxon,map$InfectionStatus=="Healthy" & map$MothChild=="Mother"]/10000,
          CPMgenus[taxon,map$InfectionStatus=="RespiratoryIllness" & map$MothChild=="Mother"]/10000,
          main=paste(taxon,"abundance in mothers"),names=c("Healthy","LRTI"))
}

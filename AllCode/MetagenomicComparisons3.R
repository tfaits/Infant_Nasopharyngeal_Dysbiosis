args <- commandArgs(trailingOnly=TRUE)
#Usage: Rscript MetagenomicComparisons2.R <FeatCountTable.txt> <OutFileBase>
FeatureCounts <- read.table(args[1],row.names=1,header=TRUE,comment.char="")
IDs <- rownames(FeatureCounts)
FeatureCounts <- data.frame(lapply(FeatureCounts,function(x){as.numeric(as.character(x))}))
rownames(FeatureCounts) <- IDs
FullConvert <- readRDS("PathoConvert2.RDS")
FullConvert$kingdom <- as.factor(FullConvert$kingdom)
FullConvert$phylum<- as.factor(FullConvert$phylum)
FullConvert$class <- as.factor(FullConvert$class)
FullConvert$order <- as.factor(FullConvert$order)
FullConvert$family <- as.factor(FullConvert$family)
FullConvert$genus <- as.factor(FullConvert$genus)
FullConvert$species <- as.factor(paste(FullConvert$genus,FullConvert$species,sep="_"))
rownames(FullConvert) <- FullConvert$IDnum

convertToTaxonLevel <- function(countMatrix, taxlevel){
   outTable <- matrix(data=NA,nrow=length(levels(FullConvert[,taxlevel])),ncol=dim(countMatrix)[2])
   for(i in 1:dim(outTable)[1]){
     outTable[i,] <- apply(countMatrix[which(FullConvert[rownames(countMatrix),taxlevel] == levels(FullConvert[,taxlevel])[i]),],2,sum)
   }
   rownames(outTable) <- levels(FullConvert[,taxlevel])
   colnames(outTable) <- colnames(countMatrix)
   return(outTable)
}

FC_p <- convertToTaxonLevel(FeatureCounts,"phylum")
saveRDS(FC_p,file=paste(args[2],"_p.RDS",sep=""))
FC_c <- convertToTaxonLevel(FeatureCounts,"class")
saveRDS(FC_c,file=paste(args[2],"_c.RDS",sep=""))
FC_o <- convertToTaxonLevel(FeatureCounts,"order")
saveRDS(FC_o,file=paste(args[2],"_o.RDS",sep=""))
FC_f <- convertToTaxonLevel(FeatureCounts,"family")
saveRDS(FC_f,file=paste(args[2],"_f.RDS",sep=""))
FC_g <- convertToTaxonLevel(FeatureCounts,"genus")
saveRDS(FC_g,file=paste(args[2],"_g.RDS",sep=""))
FC_s <- convertToTaxonLevel(FeatureCounts,"species")
saveRDS(FC_s,file=paste(args[2],"_s.RDS",sep="")) 

require(SpiecEasi)
require(igraph)

setwd("/Users/labadmin2/Documents/infant_gill/ProperLabels")
# maxCors <- c()
# for(i in 1:990){
#   fakeDat <- matrix(nrow=262,ncol=413,
#                     data=sample(genus[apply(genus[,map$Sample[map$Age<1000]],1,sum)>0,map$Sample[map$Age<1000]],
#                                 size=262*413, replace=TRUE))
#   tmp <- sparcc(fakeDat/apply(fakeDat,1,sum)*1000000, iter=20, inner_iter=10, th=0.1)
#   tmp2 <- unlist(tmp$Cor)
#   maxCors <- c(maxCors,max(abs(tmp2[tmp2<1])))
# }
maxCors <- readRDS("maxCors.RDS")


myCor <- list()
tmp <- sparcc(t(CPMgenus[apply(CPMgenus[,map$Sample[map$Age<1000]],1,sum)>0,map$Sample[map$Age<1000]]),
       iter = 20, inner_iter = 10, th = 0.1)
myCor$Gen <-tmp$Cor
for(i in 1:nrow(myCor$Gen)){
  myCor$Gen[i,i] <- 0
}
colnames(myCor$Gen) <- rownames(myCor$Gen) <- rownames(CPMgenus[apply(CPMgenus[,map$Sample[map$Age<1000]],1,sum)>0,map$Sample[map$Age<1000]])

upperGenCor <- myCor$Gen
for(i in 2:nrow(upperGenCor)){
  for(j in 1:i){
    upperGenCor[i,j] <- 0
  }
}

genCor <- melt(upperGenCor)
genCor$pVal <- unlist(lapply(genCor$value,function(x){sum(maxCors >= abs(x))/1000}))
genCor$qVal <- 1-p.adjust(genCor$pVal,method="fdr")
colnames(genCor) <- c("from","to","value","pVal","weight")
genCor$PosNeg <- "Positive"
for(i in 1:nrow(genCor)){
  if(genCor$value[i] < 0){
    genCor$PosNeg[i] <- "Negative"
  }
}
taxTable <- readRDS("PAthoConvert2.RDS")

for(i in 1:nrow(genCor)){
  genCor$Phylum[i] <- taxTable[taxTable$genus==]
}


genGraph <- graph_from_data_frame(genCor[genCor$weight > 0,], directed=FALSE)
edgeColors <- c(Positive = "darkblue", Negative = "darkred")
E(genGraph)$width <- 4
E(genGraph)$color = edgeColors[E(genGraph)$PosNeg]
plot(genGraph)

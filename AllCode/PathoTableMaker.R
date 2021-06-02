#install.packages("rentrez")
require(rentrez)
require(XML)

PathoConvert <- data.frame(IDnum="ti|1000562|org|Streptococcus_phocae_subsp._salmonis|accession|NZ_JSAP01000001.1",kingdom="Bacteria",phylum="Firmicutes",class="Bacilli",
                           order="Lactobacillales",family="Streptococcaceae",
                           genus="Streptococcus",species="phocae", stringsAsFactors=FALSE)
pullTids <- function(inTable, PathoConvert){
  for(j in 1:dim(inTable)[1]){
    if (!(rownames(inTable)[j] %in% PathoConvert$IDnum)){
      tid <- strsplit(rownames(inTable)[j],split="|",fixed=TRUE)[[1]][2]
      SpecName <- strsplit(strsplit(rownames(inTable)[j],split="|",fixed=TRUE)[[1]][4],split="_",fixed=TRUE)[[1]][2]
      r_fetch <- entrez_fetch(db="taxonomy",id=tid,rettype="xml")
      dat <- xmlToList(r_fetch)
      dat <- dat$Taxon$LineageEx
      kingdom <- "undefined"
      phylum <- "undefined"
      class <- "undefined"
      order <- "undefined"
      family <- "undefined"
      genus <- "undefined"
      species <- "undefined"
      for(i in 1:length(dat)){
        if(any(grep("kingdom",dat[[i]]$Rank))){
          kingdom <- dat[[i]]$ScientificName
        }
        if(any(grep("phylum",dat[[i]]$Rank))){
          phylum <- dat[[i]]$ScientificName
        }
        if(any(grep("class",dat[[i]]$Rank))){
          class <- dat[[i]]$ScientificName
        }
        if(any(grep("order",dat[[i]]$Rank))){
          order <- dat[[i]]$ScientificName
        }
        if(any(grep("family",dat[[i]]$Rank))){
          family <- dat[[i]]$ScientificName
        }
        if(any(grep("genus",dat[[i]]$Rank))){
          genus <- dat[[i]]$ScientificName
        }
        if(any(grep("species",dat[[i]]$Rank))){
          species <- SpecName
        } else {
          species <- SpecName
        }
      }
      PathoConvert <- rbind(PathoConvert,c(rownames(inTable)[j],kingdom,phylum,class,order,family,genus,species))
    }    
  }
  rownames(PathoConvert) <- PathoConvert$IDnum
  return(PathoConvert)
}
tmp <- read.table("../PathoscopeRuns/PathoScopeTable.txt")
ConverTable <- pullTids(tmp[1:10,],PathoConvert)
ConverTable <- pullTids(tmp[11:100,],ConverTable)
for(i in 1:28){
  ConverTable <- pullTids(tmp[(i*100 + 1):((i+1)*100),],ConverTable)
}
ConverTable <- pullTids(tmp[2900:2980,],ConverTable)
saveRDS(ConverTable,file="PathoConvert2.RDS")

require(reshape)
library(hashmap)
library(xlsx)

# This script reads the DAVID GO output chart for up and down regulated genes and add a new column
# The new column has the offician gene symbols mapped from ensembl gene IDs

path <- "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/DAVID_enrichment"
setwd(path)

files <- c("5PDF_vs_6PDF_up_DAVID_chart.txt", "5PDF_vs_6PDF_down_DAVID_chart.txt")
output <- "5PDF_vs_6PDF_GO.xlsx"

ensembleIdMap <- read.table(file = "E:/Chris_UM/Analysis/CoreData/15_ZhuBo_RNASeq2/DAVID_enrichment/ensembl_to_geneSymbol.txt", sep = "\t", header = T)

idMap <- hashmap(keys = as.vector(ensembleIdMap$ensembl_gene_id), values = as.vector(ensembleIdMap$gene_symbol))


for(i in files){
  go <- read.table(file = i, sep = "\t", header = T)
  
  sheet <- gsub(".txt", "", i)
  
  #conver the Genes columns to character
  go$Genes <- as.character(go$Genes)
  
  #map the ensembl gene id to official gene symbol
  go$Gene_symbol <- mapply(function(x){paste(idMap[[unlist(strsplit(as.character(x), split = ", "))]], collapse = ", ")}, go$Genes)
  
  write.xlsx(x = go, sheetName = sheet, file = output, col.names = T, row.names = F, append = T)
}



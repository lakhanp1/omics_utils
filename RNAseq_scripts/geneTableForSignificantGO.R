
library(dplyr)
library(tibble)
library(data.table)
require(XLConnect)


rm(list = ls())

path = "E:/Chris_UM/Analysis/CoreData/15_ZhuBo_RNASeq2/PG_PV_diff"
setwd(path)


##Prepare input data
compare <- c("PG", "PV")

goFile = "go_for_boxplot.tab"
diffFile <- "PG_vs_PV_diffData.txt"
normCountFile <- "PG_vs_PV_normCounts.tab"

outFile = "genes_specific_to_GO_boxplot.tab"


goData = fread(goFile, sep = "\t", header = T, stringsAsFactors = F, data.table = F)
diffData <- fread(diffFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)
normCounts <- fread(normCountFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)



goTermGenes = goData %>% tidyr::unnest(geneID = strsplit(Genes, ", ")) %>%
  dplyr::select(-Genes) %>%
  left_join(y = diffData, by = c("geneID")) %>%
  left_join(y = normCounts, by = c("geneID"))

fwrite(x = goTermGenes, file = outFile, sep = "\t", quote = F, col.names = T)


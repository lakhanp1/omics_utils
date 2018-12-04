
library(dplyr)
library(tibble)
library(data.table)
require(XLConnect)
options(java.parameters = "- Xmx4g")
xlcFreeMemory()


rm(list = ls())

path = "E:/Chris_UM/Analysis/CoreData/15_ZhuBo_RNASeq2/PG_PV_diff"
setwd(path)


##Prepare input data
compare <- c("PG", "PV")

geneList = "specificGeneList.tab"
diffFile <- "PG_vs_PV_diffData.txt"
normCountFile <- "PG_vs_PV_normCounts.tab"

excelOut = "PG_vs_PV_enrichment_final.xlsx"

exc = loadWorkbook(excelOut , create = FALSE)
xlcFreeMemory()


genes = fread(geneList, sep = "\t", header = T, stringsAsFactors = F, data.table = F)
diffData <- fread(diffFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)
normCounts <- fread(normCountFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)

wrkSheet = "specific_genes"
createSheet(exc, name = wrkSheet)
createFreezePane(exc, sheet = wrkSheet, 2, 2)


geneData = left_join(x = genes, y = diffData, by = c("geneID")) %>%
  left_join(y = normCounts, by = c("geneID"))

writeWorksheet(object = exc, data = geneData, sheet = wrkSheet, header = T)

xlcFreeMemory()
saveWorkbook(exc)

library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)
library(scales)
require(XLConnect)
options(java.parameters = "- Xmx4g")
xlcFreeMemory()


rm(list = ls())

source(file = "E:/Chris_UM/Codes/GO_enrichment/topGO_functions.R")

path = "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN"
setwd(path)

## This script run GO BP enrichment analysis using topGO tool for all the clusters

##################################################################################

## get the data
TF_profile = "An_laeA_20h_HA"
polII_sample = "An_untagged_20h_polII"
name = "An_laeA_20h_HA"

mapFile = "E:/Chris_UM/Database/A_Nidulans/geneid2go.ANidulans.20171004.map"

TF_dataPath = paste0("TF_data/", TF_profile, collapse = "")

clusterFile = paste0(TF_dataPath, "/", name, "_allGenes_clusters.tab", collapse = "")

topGO_path = paste0(TF_dataPath, "/", "topGO", collapse = "")

outPreTopgo_all = paste0(topGO_path, "/", name, "_allGenes", collapse = "")
outPreTopgo_expressed = paste0(topGO_path, "/", name, "_expressedGenes", collapse = "")
outPreTopgo_sm = paste0(topGO_path, "/", name, "_SM_genes", collapse = "")
outPreTopgo_peaks = paste0(topGO_path, "/", name, "_peaksGenes", collapse = "")
outPreTopgo_pkExp = paste0(topGO_path, "/", name, "_pkExpGenes", collapse = "")

polII_expId = paste("is_expressed(", polII_sample, ")", sep = "")
hasPeakCol = paste("hasPeak(", TF_profile, ")", sep = "")


clusterData = fread(file = clusterFile, sep = "\t", 
                    header = T, stringsAsFactors = F, na.strings = "NA", data.table = F)

if(!dir.exists(topGO_path)){
  dir.create(topGO_path)
}


excelOut = paste0(topGO_path, "/", name, "_topGO.xlsx", collapse = "")
unlink(excelOut, recursive = FALSE, force = FALSE)
exc = loadWorkbook(excelOut , create = TRUE)


# ## for testing purpose with one cluster
# tmpData = expressedDf %>% filter(cluster == "Cluster_7")
# 
# ## get GO enrichment table
# goData = get_topGO_enrichment(goMapFile = mapFile, genesOfInterest = tmpData$gene)
# 
# topGoScatter = topGO_scatterPlot(df = goData, title = "test plot")
# 
# imageWd = (min(max(nchar(as.character(goData$Term))), 80) * 30) * 1.5
# imageHt = max(nrow(goData) * 90, 1500)
# # imageRes = max((imageWd * imageHt / 45000), 250)
# imageRes = max(min(imageWd , imageHt) / 12, 200)
# 
# png(filename = "topGO_scatter.png", width = imageWd, height = imageHt, res = imageRes)
# print(topGoScatter)
# dev.off()



## function to generate the plot. to be called inside do()
topGO_and_plot_asDf = function(gn, tt, cl, pfx){
  tt = paste(tt, cl, sep = "\n")
  
  goData = get_topGO_enrichment(goMapFile = mapFile, genesOfInterest = gn)
  
  if(nrow(goData) == 0){
    return(data.frame(height = NA, width = NA, title = NA, res = NA, png = NA, stringsAsFactors = F))
  }
  
  topGoScatter = topGO_scatterPlot(df = goData, title = tt)
  
  ht = max(nrow(goData) * 80, 1500)
  wd = (min(max(nchar(as.character(goData$Term))), 80) * 30) * 1.5
  rs = max(min(wd, ht) / 12, 200)
  
  pngFile = paste(pfx, "_", cl, "_topGO.png", sep = "")
  
  png(filename = pngFile, width = wd, height = ht, res = rs)
  print(topGoScatter)
  dev.off()
  
  return(data.frame(height = ht, 
                    width = wd, 
                    res = rs, 
                    count = nrow(goData), 
                    title = tt, 
                    png = pngFile,
                    stringsAsFactors = F))
}


##################################################################################
## run topGO: all genes

allTitle = paste("GO enrichment using topGO for all genes \n TF:", TF_profile, "and polII:", polII_sample, sep = " ")

## GO enrichment table for all the clusters
goEnrichment_all = clusterData %>% 
  group_by(cluster) %>%
  do(get_topGO_enrichment(goMapFile = mapFile, genesOfInterest = .$gene))


fwrite(goEnrichment_all, file = paste(outPreTopgo_all, "_topGO.tab"),
       sep = "\t", row.names = F, col.names = T, quote = F)


## write data to Excel
xlcFreeMemory()
wrkSheet = "allGenes"
createSheet(exc, name = wrkSheet)
createFreezePane(exc, sheet = wrkSheet, 2, 2)
setMissingValue(object = exc, value = "NA")
writeWorksheet(object = exc, data = goEnrichment_all, sheet = wrkSheet, header = T)
setAutoFilter(object = exc, sheet = wrkSheet, reference = aref(topLeft = "A1", dimension = dim(goEnrichment_all)))
xlcFreeMemory()



## generate scatter plots for each
topGOPlots = clusterData %>% 
  group_by(cluster) %>%
  do(topGO_and_plot_asDf(gn = .$gene, 
                         tt = allTitle, 
                         cl = unique(.$cluster),
                         pfx = outPreTopgo_all
                         )
     )

fwrite(topGOPlots, file = paste(outPreTopgo_all, "_topGoStats.tab"),
       sep = "\t", row.names = F, col.names = T, quote = F)




##################################################################################
## run topGO: top 10% expressed genes
expTitle = paste("GO enrichment using topGO for polII expressed genes \n TF:", TF_profile, "and polII:", polII_sample, sep = " ")

expressedDf = clusterData %>% filter(UQ(as.name(polII_expId)) == "TRUE")

## GO enrichment table for all the clusters
goEnrichment_exp = expressedDf %>% 
  group_by(cluster) %>%
  do(get_topGO_enrichment(goMapFile = mapFile, genesOfInterest = .$gene))


fwrite(goEnrichment_exp, file = paste(outPreTopgo_expressed, "_topGO.tab"),
       sep = "\t", row.names = F, col.names = T, quote = F)


## write data to Excel
xlcFreeMemory()
wrkSheet = "expressedGenes"
createSheet(exc, name = wrkSheet)
createFreezePane(exc, sheet = wrkSheet, 2, 2)
setMissingValue(object = exc, value = "NA")
writeWorksheet(object = exc, data = goEnrichment_exp, sheet = wrkSheet, header = T)
setAutoFilter(object = exc, sheet = wrkSheet, reference = aref(topLeft = "A1", dimension = dim(goEnrichment_exp)))
xlcFreeMemory()


## generate scatter plots for each
topGOPlots_exp = expressedDf %>% 
  group_by(cluster) %>%
  do(topGO_and_plot_asDf(gn = .$gene, 
                         tt = expTitle, 
                         cl = unique(.$cluster),
                         pfx = outPreTopgo_expressed
  )
  )

fwrite(topGOPlots_exp, file = paste(outPreTopgo_expressed, "_topGoStats.tab"),
       sep = "\t", row.names = F, col.names = T, quote = F)




##################################################################################
## run topGO: genes for which peak was called by macs2

peakTitle = paste("GO enrichment using topGO for genes bound by TF \n TF:", TF_profile, "and polII:", polII_sample, sep = " ")

peaksDf = clusterData %>% filter(UQ(as.name(hasPeakCol)) == "TRUE")

## GO enrichment table for all the clusters
goEnrichment_peak = peaksDf %>% 
  group_by(cluster) %>%
  do(get_topGO_enrichment(goMapFile = mapFile, genesOfInterest = .$gene))


fwrite(goEnrichment_peak, file = paste(outPreTopgo_peaks, "_topGO.tab"),
       sep = "\t", row.names = F, col.names = T, quote = F)


## write data to Excel
xlcFreeMemory()
wrkSheet = "tfTargetGenes"
createSheet(exc, name = wrkSheet)
createFreezePane(exc, sheet = wrkSheet, 2, 2)
setMissingValue(object = exc, value = "NA")
writeWorksheet(object = exc, data = goEnrichment_peak, sheet = wrkSheet, header = T)
setAutoFilter(object = exc, sheet = wrkSheet, reference = aref(topLeft = "A1", dimension = dim(goEnrichment_peak)))
xlcFreeMemory()



## generate scatter plots for each
topGOPlots_peak = peaksDf %>% 
  group_by(cluster) %>%
  do(topGO_and_plot_asDf(gn = .$gene, 
                         tt = peakTitle, 
                         cl = unique(.$cluster),
                         pfx = outPreTopgo_peaks
  )
  )

fwrite(topGOPlots_peak, file = paste(outPreTopgo_peaks, "_topGoStats.tab"),
       sep = "\t", row.names = F, col.names = T, quote = F)





##################################################################################
## run topGO: polII expressed and TF bound genes together

pkExp_title = paste("GO enrichment using topGO for genes bound by TF and expressed in WT polII\n TF:", TF_profile, "and polII:", polII_sample, sep = " ")

## select the genes which are expressed in polII sample OR have TSS peak
pkExpdf = filter_at(.tbl = clusterData, .vars = c(polII_expId, hasPeakCol), .vars_predicate = any_vars(. == "TRUE"))

## GO enrichment table for all the clusters
goEnrichment_pkExp = pkExpdf %>% 
  group_by(cluster) %>%
  do(get_topGO_enrichment(goMapFile = mapFile, genesOfInterest = .$gene))


fwrite(goEnrichment_pkExp, file = paste(outPreTopgo_pkExp, "_topGO.tab"),
       sep = "\t", row.names = F, col.names = T, quote = F)



## write data to Excel
xlcFreeMemory()
wrkSheet = "peakExpressedGenes"
createSheet(exc, name = wrkSheet)
createFreezePane(exc, sheet = wrkSheet, 2, 2)
setMissingValue(object = exc, value = "NA")
writeWorksheet(object = exc, data = goEnrichment_pkExp, sheet = wrkSheet, header = T)
setAutoFilter(object = exc, sheet = wrkSheet, reference = aref(topLeft = "A1", dimension = dim(goEnrichment_pkExp)))
xlcFreeMemory()



## generate scatter plots for each
topGOPlots_pkExp = pkExpdf %>% 
  group_by(cluster) %>%
  do(topGO_and_plot_asDf(gn = .$gene, 
                         tt = pkExp_title, 
                         cl = unique(.$cluster),
                         pfx = outPreTopgo_pkExp
  )
  )

fwrite(topGOPlots_pkExp, file = paste(outPreTopgo_pkExp, "_topGoStats.tab"),
       sep = "\t", row.names = F, col.names = T, quote = F)







saveWorkbook(exc)




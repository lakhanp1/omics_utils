library(dplyr)
library(data.table)
library(tibble)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)
library(amap)
library(imputeTS)
library(lazyeval)
require(XLConnect)
options(java.parameters = "- Xmx4g")
xlcFreeMemory()

## This script plots the profile heatmaps for multiple samples. If the expression values are present, it
## also plots a simple heatmap using these expression values
## Original script: multiProfilePlots.R
## Modified for Fangfang

rm(list = ls())

source(file = "E:/Chris_UM/Codes/ChIP_seq/multiTF_analysis_functions.R")
source(file = "E:/Chris_UM/Codes/ChIP_seq/TF_polII_functions.R")

path = "E:/Chris_UM/Codes/ChIP_seq/profilePlotsFangFang/data"
setwd(path)


##################################################################################
## configurations
comparisonName = "testFang"
clusterNumbers = 8
generateClusters = TRUE
showExpressionHeatmap = FALSE


outDir = path

geneInfoFile = "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"
clusterOrderFile = "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/referenceData/clusterOrder.txt"
anGenes = "E:/Chris_UM/Codes/ChIP_seq/profilePlotsFangFang/AN_genesForPolII.bed"

sampleFile = "E:/Chris_UM/Codes/ChIP_seq/profilePlotsFangFang/fangfangProfileInput.txt"


##################################################################################


outPrefix_all = paste0(outDir, "/", comparisonName, "_allGenes", collapse = "")
outPrefix_select = paste0(outDir, "/", comparisonName, "_selectedGenes", collapse = "")
outPrefix_sm = paste0(outDir, "/", comparisonName, "_SM_genes", collapse = "")
outPrefix_peaks = paste0(outDir, "/", comparisonName, "_peaksGenes", collapse = "")
outPrefix_pkExp = paste0(outDir, "/", comparisonName, "_pkExpGenes", collapse = "")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}



## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
exptData = fread(file = sampleFile, sep = "\t", stringsAsFactors = F, header = T)

## genes to read
geneSet = fread(file = anGenes, header = F, select = c(4), col.names = c("gene"), stringsAsFactors = F)

##################################################################################

## read the profile matrix
# mat1 = getDeeptoolsProfileMatrix(file = exptData$matFile[1], 
#                                  signalName = exptData$Sample_ID[1],
#                                  selectGenes = geneSet)

mat1 = getMiaoProfileMatrix(file = exptData$matFile[1],
                            signalName = exptData$Sample_ID[1],
                            selectGenes = geneSet, up = 100, target = 50, down = 50)


## check the distribution in data
quantile(mat1, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)
tfColFun = colorRamp2(quantile(mat1, c(0.50, 0.99), na.rm = T), c("white", "red"))
polIIColFun = colorRamp2(quantile(mat1, c(0.01, 0.5, 0.99), na.rm = T), c("blue", "white", "red"))

col_fun = tfColFun
if(exptData$IP_tag[1] == "polII"){
  col_fun = polIIColFun
}

## get first primary profile heatmap
profile1 = primaryProfileHeatmap(profileMat = mat1, 
                                 signalName = exptData$profileName[1], 
                                 numClust = clusterNumbers, 
                                 color = col_fun,
                                 createClusters = generateClusters,
                                 clusterStorePath = exptData$clusterFile[1])


##################################################################################

## gene information annotations: cluster and TF and polII expression values
geneInfo = getGeneInfo(dataFile = geneInfoFile, orderFile = clusterOrderFile)

geneInfo = left_join(x = profile1$cluster, geneInfo, by = c("gene" = "gene"))

head(geneInfo)

hasPeakCol = "has_TSS_peak"
polII_ids = exptData$Sample_ID[which(exptData$IP_tag == "polII")]
polII_expIds = paste("is_expressed(", polII_ids, ")", sep = "")

expressionData = get_polII_expressions(exptInfo = exptData,
                                       genesDf = geneInfo)


expressionData = get_TF_binding_data(exptInfo = exptData, 
                                     genesDf = expressionData, 
                                     peakCol = hasPeakCol)

##################################################################################


multiProfiles_all = multiProfilePlots(exptInfo = exptData, 
                                      expressionData = expressionData, 
                                      genesToPlot = geneSet,
                                      clusters = profile1$cluster,
                                      clusterColor = profile1$clusterColor,
                                      matSource = "miao",
                                      matBins = c(100, 50, 50),
                                      profileColors = NULL, 
                                      expressionColor = NULL, 
                                      plotExpression = showExpressionHeatmap)



wd = (nrow(exptData) * 3000) + (length(polII_ids) * 1000 * showExpressionHeatmap)
title_all = "Transcription factor and polII binding profile"

# draw Heatmap and add the annotation name decoration
png(filename = paste(outPrefix_all, "_profiles.png", sep = ""), width=wd, height=10000, res = 800)

draw(multiProfiles_all$heatmapList, 
     main_heatmap = exptData$profileName[1],
     annotation_legend_list = list(profile1$legend),
     column_title = title_all,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(multiProfiles_all$plotGaps, "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()


##################################################################################

## plot only genes of interest

## select the genes: copy from clipboard
geneSubset = data.frame(gene = readClipboard(), stringsAsFactors = F)

subsetDf = left_join(x = geneSet, y = expressionData, by = c("gene" = "gene"))
row.names(subsetDf) = subsetDf$gene

clusterInfo = subsetDf["cluster"]


multiProfiles_select = multiProfilePlots(exptInfo = exptData, 
                                         expressionData = subsetDf, 
                                         genesToPlot = subsetDf["gene"],
                                         clusters = clusterInfo,
                                         clusterColor = profile1$clusterColor,
                                         matSource = "miao",
                                         matBins = c(100, 50, 50),
                                         profileColors = multiProfiles_all$profileColors, 
                                         expressionColor = multiProfiles_all$expressionColor,
                                         plotExpression = showExpressionHeatmap)


wd = (nrow(exptData) * 3000) + (length(polII_ids) * 1000 * showExpressionHeatmap)
title_exp = "Transcription factor and polII binding profile: polII expressed genes"

# draw Heatmap and add the annotation name decoration
png(filename = paste(outPrefix_select, "_profiles.png", sep = ""), width=wd, height=10000, res = 800)

draw(multiProfiles_select$heatmapList, 
     main_heatmap = exptData$profileName[1],
     annotation_legend_list = list(profile1$legend),
     column_title = title_exp,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(multiProfiles_select$plotGaps, "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()




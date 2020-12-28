suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))


## This script plots the profile heatmaps for multiple samples. If the expression values are present, it
## also plots a simple heatmap using these expression values

rm(list = ls())

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "geneset2_sex_asex_20h"
workDir <- here::here("analysis", "03_KERS_complex", "KERS_complex_20h", comparisonName)
outPrefix <- paste(workDir, "/", comparisonName, sep = "")

file_plotSamples <- paste(workDir, "/", "samples.txt", sep = "")
file_geneSubset <- paste(workDir, "/", "geneList.txt", sep = "")

matrixType <- "2kb_ATG_1kb"
up <- 2000
body <- 0
down <- 1000
binSize <- 10

matrixDim = c(c(up, body, down)/binSize, binSize)

showExpressionHeatmap <- FALSE

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF


## colors
colList <- list()

outPrefix_all <- paste(outPrefix, "_allGenes", sep = "")
outPrefix_geneset <- paste(outPrefix, "_genes", sep = "")


##################################################################################
sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, comment = "#"))

tempSInfo <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = sampleList$sampleId,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)

tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST")]


## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = tfIds,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)


inputData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = inputIds,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = polII_ids,
  dataPath = polII_dataPath, profileMatrixSuffix = matrixType
)

histData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = histIds,
  dataPath = hist_dataPath, profileMatrixSuffix = matrixType
)


exptData <- dplyr::bind_rows(tfData, inputData, histData, polIIData)

exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)

##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand"))

kmClust <- dplyr::left_join(
  x = suppressMessages(readr::read_tsv(file = tfData$clusterFile[1])),
  y = geneSet, by = c("geneId" = "geneId")
)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, keytype = "GID",
                                  columns = c("GENE_NAME", "DESCRIPTION")) %>% 
  dplyr::rename(geneId = GID)

# geneInfo <- dplyr::left_join(x = kmClust, y = geneDesc, by = c("geneId" = "GID"))
geneInfo <- geneDesc

head(geneInfo)


expressionData <- get_TF_binding_data(exptInfo = tfData,
                                      genesDf = geneInfo)

# expressionData <- get_polII_expressions(exptInfo = polIIData,
#                                         genesDf = expressionData)

# view(dfSummary(expressionData))

peakTargetMat <- peak_target_matrix(sampleInfo = tfData, position = "best")


anLables <- list()
# anLables[[tssPeakTypeCol]] = gsub("peakType", "TSS peak type\n", tssPeakTypeCol) %>% gsub("\\(|\\)", "", .)
# anLables[[tesPeakTypeCol]] = gsub("tesPeakType", "TES peak type\n", tesPeakTypeCol) %>% gsub("\\(|\\)", "", .)
# anLables[[isExpCol]] = txt = gsub("is_expressed", "is expressed\n", isExpCol) %>% gsub("\\(|\\)", "", .)
anLables[["is_SM_gene"]] = "SM gene"
anLables[["is_TF"]] = "Transcription Factor"
anLables[["gene_length"]] = "Gene Length"

##################################################################################


## generate profile matrix for the first time
# genesGr <- rtracklayer::import(con = file_genes, format = "bed")
# 
# geneStartGr <- GenomicRanges::resize(x = genesGr, width = 1, fix = "start")
# 
# i <- 1
# for(i in 1:nrow(exptData)){
#   bwMat <- bigwig_profile_matrix(bwFile = exptData$bwFile[i],
#                                  regions = geneStartGr,
#                                  signalName = exptData$sampleId[i],
#                                  extend = c(up, down),
#                                  targetName = "ATG",
#                                  storeLocal = TRUE,
#                                  localPath = exptData$matFile[i])
# }

## color list
matList <- import_profiles(exptInfo = dplyr::bind_rows(tfData, inputData, histData, polIIData),
                           geneList = geneInfo$geneId,
                           targetType = "point", targetName = "ATG",
                           up = matrixDim[1], target = matrixDim[2], down = matrixDim[3])


## tf colors
tfMeanProfile <- NULL
if(length(c(tfIds)) == 1){
  tfMeanProfile <- matList[[tfIds]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[tfIds])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))
tfColorList <- sapply(
  X = c(tfIds, inputIds),
  FUN = function(x){
    return(
      colorRamp2(breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
                 colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))))
  }
)

## polII colors
polIIMeanProfile <- NULL
polIIColorList <- NULL
# if(nrow(polIIData) == 1){
#   polIIMeanProfile <- matList[[polIIData$sampleId]]
# } else{
#   polIIMeanProfile <- getSignalsFromList(lt = matList[polIIData$sampleId])
# }
# quantile(polIIMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# polIIMeanColor <- colorRamp2(quantile(polIIMeanProfile, c(0.01, 0.5, 0.995), na.rm = T), c("blue", "white", "red"))
# polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})

## histone colors
histMeanProfile <- NULL
histColorList <- NULL
# if(nrow(histData) == 1){
#   histMeanProfile <- matList[[histData$sampleId]]
# } else{
#   histMeanProfile <- getSignalsFromList(lt = matList[histData$sampleId])
# }
# quantile(histMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# # histMeanColor <- colorRamp2(quantile(histMeanProfile, c(0.30, 0.995), na.rm = T), c("black", "yellow"))
# histColorList <- sapply(
#   X = histData$sampleId,
#   FUN = function(x){
#     return(colorRamp2(breaks = quantile(histMeanProfile, c(0.30, 0.995), na.rm = T),
#                       colors =  unlist(strsplit(x = exptDataList[[x]]$color, split = ","))))
#   }
# )

colorList <- unlist(list(tfColorList, polIIColorList, histColorList))
ylimList <- list()
# ylimList <- sapply(c(polII_ids, histIds), function(x){return(0.996)}, simplify = FALSE)
ylimList <- append(x = ylimList,
                   values = sapply(c(tfIds, inputIds), function(x){return(c(0, 25))}, simplify = FALSE))



##################################################################################

# 
# multiProfiles_all <- multi_profile_plots(exptInfo = exptData,
#                                          genesToPlot = geneInfo$geneId,
#                                          targetType = "point",
#                                          targetName = "ATG",
#                                          matBins = matrixDim,
#                                          clusters = expressionData,
#                                          profileColors = colorList,
#                                          column_title_gp = gpar(fontsize = 12),
#                                          ylimFraction = ylimList)
# 
# 
# ## gene length annotation
# anGl <- gene_length_heatmap_annotation(bedFile = file_genes, genes = expressionData$geneId)
# 
# 
# htlist_allGenes <- anGl$an +
#   multiProfiles_all$heatmapList
# # tssTesPlotList
# 
# # ## row order as decreasing polII signal
# # if( all(rownames(htlist_allGenes@ht_list[[ polIIData$profileName[1] ]]@matrix) == expressionData$gene) ){
# #   
# #   rowOrd <- order(expressionData[[polIICols$exp[1]]], decreasing = TRUE)
# # }
# 
# pdfWd <- 2 + 
#   (length(multiProfiles_all$heatmapList@ht_list) * 2) +
#   (length(polII_ids) * 0.25 * showExpressionHeatmap)
# 
# title_all <- paste(comparisonName, ": all genes", collapse = "")
# 
# # draw Heatmap and add the annotation name decoration
# # png(filename = paste(outPrefix_all, ".profiles.png", sep = ""), width=wd, height=3500, res = 250)
# pdf(file = paste(outPrefix_all, ".profiles.pdf", sep = ""), width = pdfWd, height = 13)
# 
# draw(htlist_allGenes,
#      main_heatmap = exptData$profileName[1],
#      # annotation_legend_list = list(profile1$legend),
#      column_title = title_all,
#      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
#      row_sub_title_side = "left",
#      heatmap_legend_side = "bottom",
#      gap = unit(7, "mm"),
#      # row_order = rowOrd,
#      padding = unit(rep(0.5, times = 4), "cm")
# )
# 
# dev.off()

##################################################################################
# plot profiles for genes of interest

geneSubset <- suppressMessages(
  readr::read_tsv(file = file_geneSubset)
)

geneSubset <- dplyr::left_join(x = geneSubset, y = peakTargetMat, by = "geneId") %>% 
  dplyr::left_join(y = geneInfo, by = "geneId")


multiProfiles_geneset <- multi_profile_plots(
  exptInfo = exptData,
  genesToPlot = geneSubset$geneId,
  targetType = "point",
  targetName = "ATG",
  matBins = matrixDim,
  clusters = NULL,
  showAnnotation = FALSE,
  profileColors = colorList,
  column_title_gp = gpar(fontsize = 12),
  # row_order = geneSubset$geneId,
  # show_row_names = TRUE,
  # row_labels = geneSubset$geneName,
  ylimFraction = ylimList
)


pdfWd <- 2 + 
  (length(multiProfiles_geneset$heatmapList@ht_list) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap)


## gene length annotation
anGl_geneset <- gene_length_heatmap_annotation(
  bedFile = file_genes,
  genes = geneSubset$geneId,
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb")),
  pointSize = unit(4, "mm"))


geneset_htlist <- anGl_geneset$an +
  multiProfiles_geneset$heatmapList



# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_geneset = paste(comparisonName, ": genes of interest", collapse = "")

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix_geneset, ".profiles.pdf", sep = ""), width = pdfWd, height = 12)

geneset_htlist <- draw(geneset_htlist,
                       main_heatmap = exptData$profileName[1],
                       # annotation_legend_list = list(profile1$legend),
                       split = geneSubset$group,
                       column_title = title_geneset,
                       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                       row_sub_title_side = "left",
                       heatmap_legend_side = "bottom",
                       gap = unit(7, "mm"),
                       padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()


rowOrderDf <- row_order(geneset_htlist) %>% 
  purrr::map_dfr(.f = ~ tibble(rank = ., geneId = geneSubset$geneId[.]))

ordered_data <- dplyr::left_join(x = rowOrderDf, y = geneSubset, by = "geneId")

readr::write_tsv(x = ordered_data, file = paste(outPrefix_geneset, ".data.tab", sep = ""))

##################################################################################












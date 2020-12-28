suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))


## This script plots the profile heatmaps for multiple samples. If the expression values are present, it
## also plots a simple heatmap using these expression values

rm(list = ls())

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "kdmB_del_sntB_del_20h_combined"
outPrefix <- here::here("analysis", "04_KS_del_binding", comparisonName, comparisonName)

file_plotSamples <- here::here("analysis", "04_KS_del_binding", comparisonName, "samples.txt")

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

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
outPrefix_expressed <- paste(outPrefix, "_expressedGenes", sep = "")
outPrefix_sm <- paste(outPrefix, "_SM_genes", sep = "")
outPrefix_peaks <- paste(outPrefix, "_peaksGenes", sep = "")
outPrefix_pkExp <- paste(outPrefix, "_pkExpGenes", sep = "")

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


# histTssData <- get_sample_information(
#   exptInfoFile = file_exptInfo, samples = histIds,
#   dataPath = hist_dataPath, profileMatrixSuffix = tssMatType,
#   profileType = "TSS_profile"
# )
# 
# histTesData <- get_sample_information(
#   exptInfoFile = file_exptInfo, samples = histIds,
#   dataPath = hist_dataPath, profileMatrixSuffix = tesMatType,
#   profileType = "TES_profile"
# )


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
                                  columns = c("GENE_NAME", "DESCRIPTION"))

smGenes <- AnnotationDbi::select(
  x = orgDb, keys = keys(x = orgDb, keytype = "SM_ID"), keytype = "SM_ID",
  columns = c("GID")
) %>% 
  dplyr::distinct(GID) %>% 
  dplyr::mutate(is_SM_gene = TRUE)

geneInfo <- dplyr::left_join(x = kmClust, y = geneDesc, by = c("geneId" = "GID")) %>% 
  dplyr::left_join(y = smGenes, by = c("geneId" = "GID"))

head(geneInfo)

expressionData <- get_TF_binding_data(exptInfo = tfData,
                                      genesDf = geneInfo)

# expressionData <- get_polII_expressions(exptInfo = polIIData,
#                                         genesDf = expressionData)

glimpse(expressionData)


anLables <- list()
# anLables[[tssPeakTypeCol]] = gsub("peakType", "TSS peak type\n", tssPeakTypeCol) %>% gsub("\\(|\\)", "", .)
# anLables[[tesPeakTypeCol]] = gsub("tesPeakType", "TES peak type\n", tesPeakTypeCol) %>% gsub("\\(|\\)", "", .)
# anLables[[isExpCol]] = txt = gsub("is_expressed", "is expressed\n", isExpCol) %>% gsub("\\(|\\)", "", .)
anLables[["is_SM_gene"]] = "SM gene"
anLables[["is_TF"]] = "Transcription Factor"
anLables[["gene_length"]] = "Gene Length"

##################################################################################
## color list
matList <- import_profiles(exptInfo = dplyr::bind_rows(tfData, inputData, histData, polIIData),
                           geneList = expressionData$geneId,
                           source = matrixType,
                           up = matrixDim[1], target = matrixDim[2], down = matrixDim[3])

# matListTss <- import_profiles(exptInfo = dplyr::bind_rows(tfTssData, histTssData),
#                                   geneList = geneInfo$gene,
#                                   source = "normalizedmatrix",
#                                   up = tssMatDim[1], down = tssMatDim[3],
#                                   targetType = "TSS")
# 
# matListTes <- import_profiles(exptInfo = dplyr::bind_rows(tfTesData, histTesData),
#                                   geneList = geneInfo$gene,
#                                   source = "normalizedmatrix",
#                                   up = tesMatDim[1], down = tesMatDim[3],
#                                   targetType = "TES")




# matrix_list_color <- function(df, matrices, colorType){
#   if(length(matrices) > 1){
#     meanProfile <- getSignalsFromList(lt = matrices[df$sampleId])
#   } else if(length(matrices) == 1){
#     meanProfile = matrices[[1]]
#   }
#   
#   profileColor <- NULL
#   
#   quantile(meanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
#   if(tolower(colorType) == "tf"){
# 
#     profileColor <- colorRamp2(quantile(meanProfile, c(0.05, 0.50, 0.995), na.rm = T),
#                                c("#4d9221", "#f7f7f7", "#c51b7d"))
# 
#   } else if(tolower(colorType) == "polii"){
# 
#     profileColor <- colorRamp2(quantile(meanProfile, c(0.01, 0.5, 0.995), na.rm = T),
#                                c("blue", "white", "red"))
# 
#   } else if(tolower(colorType) == "hist"){
#     profileColor <- colorRamp2(quantile(meanProfile, c(0.20, 0.995), na.rm = T),
#                                c("black", "yellow"))
#   }
# 
#   colorList <- sapply(X = df$sampleId, FUN = function(x){return(profileColor)})
# 
#   return(colorList)
# }

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
if(nrow(histData) == 1){
  histMeanProfile <- matList[[histData$sampleId]]
} else{
  histMeanProfile <- getSignalsFromList(lt = matList[histData$sampleId])
}
quantile(histMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# histMeanColor <- colorRamp2(quantile(histMeanProfile, c(0.30, 0.995), na.rm = T), c("black", "yellow"))
histColorList <- sapply(
  X = histData$sampleId,
  FUN = function(x){
    return(colorRamp2(breaks = quantile(histMeanProfile, c(0.30, 0.995), na.rm = T),
                      colors =  unlist(strsplit(x = exptDataList[[x]]$color, split = ","))))
  }
)

colorList <- unlist(list(tfColorList, polIIColorList, histColorList))
ylimList <- list()
# ylimList <- sapply(c(polII_ids, histIds), function(x){return(0.996)}, simplify = FALSE)
ylimList <- append(x = ylimList,
                   values = sapply(c(tfIds, inputIds), function(x){return(c(0, 25))}, simplify = FALSE))



##################################################################################

# expressionData$group <- dplyr::group_by_at(expressionData, .vars = unname(polIICols$is_expressed)) %>% 
#   dplyr::group_indices()
# 
# newClusters <- dplyr::select(expressionData, gene, group) %>% 
#   dplyr::rename(cluster = group)

multiProfiles_all <- multi_profile_plots(
  exptInfo = exptData,
  expressionData = expressionData,
  genesToPlot = expressionData$geneId,
  matSource = matrixType,
  matBins = matrixDim,
  clusters = expressionData,
  clusterColor = NULL,
  # clustOrd = rev(unique(newClusters$cluster)),
  profileColors = colorList,
  expressionColor = NULL,
  plotExpression = showExpressionHeatmap,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList
)

# tssProfiles_all <- multi_profile_plots(exptInfo = dplyr::bind_rows(tfTssData, histTssData),
#                                        expressionData = expressionData,
#                                        genesToPlot = geneInfo$gene,
#                                        targetType = "TSS",
#                                        matSource = "normalizedmatrix",
#                                        matBins = tssMatDim,
#                                        clusters = geneInfo,
#                                        clusterColor = multiProfiles_all$profileHeatmaps[[1]]$clusterColor,
#                                        drawClusterAn = FALSE,
#                                        profileColors = colorList,
#                                        column_title_gp = gpar(fontsize = 12),
#                                        ylimFraction = ylimList)
# 
# tesProfiles_all <- multi_profile_plots(exptInfo = dplyr::bind_rows(tfTesData, histTesData),
#                                        expressionData = expressionData,
#                                        genesToPlot = geneInfo$gene,
#                                        targetType = "TES",
#                                        matSource = "normalizedmatrix",
#                                        matBins = tesMatDim,
#                                        clusters = geneInfo,
#                                        clusterColor = multiProfiles_all$profileHeatmaps[[1]]$clusterColor,
#                                        drawClusterAn = FALSE,
#                                        profileColors = colorList,
#                                        column_title_gp = gpar(fontsize = 12),
#                                        ylimFraction = ylimList)
# 
# 
# tssTesPairs <- purrr::map2(.x = tssProfiles_all$profileHeatmaps,
#                            .y = tesProfiles_all$profileHeatmaps,
#                            .f =  function(x, y) x$heatmap + y$heatmap)
# 
# tssTesPlotList <- NULL
# 
# for (p in names(tssTesPairs)) {
#   tssTesPlotList <- tssTesPlotList + tssTesPairs[[p]]
# }

## gene length annotation
anGl <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = expressionData$geneId,
  # pointSize = unit(4, "mm"),
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb"))
)


htlist_allGenes <- anGl$an +
  multiProfiles_all$heatmapList
# tssTesPlotList

# ## row order as decreasing polII signal
# if( all(rownames(htlist_allGenes@ht_list[[ polIIData$profileName[1] ]]@matrix) == expressionData$gene) ){
#   
#   rowOrd <- order(expressionData[[polIICols$exp[1]]], decreasing = TRUE)
# }

pdfWd <- 2 + 
  (length(multiProfiles_all$heatmapList@ht_list) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap)

title_all <- paste(comparisonName, ": all genes", collapse = "")

# png(filename = paste(outPrefix_all, "_profiles.png", sep = ""), width=wd, height=3500, res = 250)
pdf(file = paste(outPrefix_all, "_profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(htlist_allGenes,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_all,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(7, "mm"),
     # row_order = rowOrd,
     padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()

##################################################################################
# plot genes which has TF peak in any of the TF samples
hasPeakDf <- filter_at(
  .tbl = expressionData,
  .vars = unname(tfCols$hasPeak[ sampleList$sampleId[sampleList$usePeaks] ]),
  .vars_predicate = any_vars(. == "TRUE")
)

## plot genes which has TF peak in a specific TF sample
# hasPeakDf <- filter_at(
#   .tbl = expressionData,
#   .vars = unname(tfCols$hasPeak[1]),
#   .vars_predicate = any_vars(. == "TRUE")
# )


multiProfiles_peak <- multi_profile_plots(
  exptInfo = exptData,
  expressionData = hasPeakDf,
  genesToPlot = hasPeakDf$geneId,
  matSource = matrixType,
  matBins = matrixDim,
  clusters = hasPeakDf,
  clusterColor = multiProfiles_all$profileHeatmaps[[1]]$clusterColor,
  profileColors = colorList,
  expressionColor = multiProfiles_all$expressionColor,
  column_title_gp = gpar(fontsize = 12),
  plotExpression = showExpressionHeatmap,
  ylimFraction = ylimList
)


## gene length annotation
anGl_peaks <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = hasPeakDf$geneId,
  # pointSize = unit(4, "mm"),
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb"))
)


peaks_htlist <- anGl_peaks$an +
  multiProfiles_peak$heatmapList


## make sure that the order of genes in the heatmap list and in the dataframe is same
if(all(rownames(peaks_htlist@ht_list[[ tfData$profileName[1] ]]@matrix) == hasPeakDf$geneId)){
  
  rowOrd_peaks <- order(hasPeakDf[[ tfCols$peakDist[[1]] ]], decreasing = TRUE)
  
}


# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_peak= paste(comparisonName, ": genes with binding signal (macs2 peak targets)", collapse = "")


pdf(file = paste(outPrefix_peaks, "_profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(peaks_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(7, "mm"),
     # row_order = rowOrd_peaks,
     padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()


## draw ungrouped profile Heatmap and add the annotation name decoration
# 
# peaks_htlist2 <- anGl_peaks$an
# i <- 1
# for(i in 1:nrow(exptData)){
#   peaks_htlist2 <- peaks_htlist2 + multiProfiles_peak$profileHeatmaps[[exptData$sampleId[i]]]$heatmap
# }
# 
# 
# pdf(file = paste(outPrefix_peaks, "_profiles_ungrouped.pdf", sep = ""), width = pdfWd, height = 13)
# 
# draw(peaks_htlist2,
#      main_heatmap = exptData$profileName[1],
#      # annotation_legend_list = list(profile1$legend),
#      column_title = title_peak,
#      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
#      row_sub_title_side = "left",
#      heatmap_legend_side = "bottom",
#      gap = unit(7, "mm"),
#      split = rep(1, nrow(hasPeakDf)),
#      # row_order = rowOrd_peaks,
#      padding = unit(rep(0.5, times = 4), "cm")
# )
# 
# dev.off()


##################################################################################
## plot polII expressed and TF bound genes together

## select the genes which are expressed in polII sample OR have TSS peak
peakExpDf <- filter_at(.tbl = expressionData,
                       .vars = unname(c(tfCols$hasPeak, polIICols$is_expressed)),
                       .vars_predicate = any_vars(. == "TRUE"))

dplyr::group_by_at(peakExpDf, .vars = unname(c(tfCols$hasPeak))) %>%
  dplyr::summarise(n= n())

# peakExpDf$group <- dplyr::group_by_at(peakExpDf, .vars = unname(c(tfCols$hasPeak, polIICols$is_expressed))) %>%
#   dplyr::group_indices()
# 
# newClusters <- dplyr::select(peakExpDf, gene, group) %>%
#   dplyr::rename(cluster = group)

multiProfiles_pkExp <- multi_profile_plots(
  exptInfo = exptData,
  expressionData = peakExpDf,
  genesToPlot = peakExpDf$geneId,
  matSource = matrixType,
  matBins = matrixDim,
  clusters = NULL,
  # clusterColor = multiProfiles_all$profileHeatmaps[[1]]$clusterColor,
  profileColors = colorList,
  expressionColor = multiProfiles_all$expressionColor,
  plotExpression = showExpressionHeatmap,
  ylimFraction = ylimList
)


## heatmap of binary assignment of samples to different group
htMat <- dplyr::select(peakExpDf, gene, starts_with("hasPeak")) %>% 
  dplyr::mutate_if(.predicate = is.logical, .funs = as.character) %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("gene")

## column name as annotation for Heatmap
colNameAnn <- HeatmapAnnotation(
  colName = anno_text(
    x = unname(tfCols$hasPeak),
    rot = 90, just = "left",
    offset = unit(1, "mm"),
    gp = gpar(fontsize = 10)
  )
)

grpHt_pkExp <- Heatmap(
  htMat,
  col = c("TRUE" = "black", "FALSE" = "white"),
  heatmap_legend_param = list(title = "Peak detected"),
  # column_names_side = "top",
  show_column_names = FALSE,
  top_annotation = colNameAnn,
  cluster_columns = FALSE, cluster_rows = FALSE,
  width = unit(3, "cm"),
  show_row_names = FALSE
)

## gene length annotation
anGl_pkExp <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = peakExpDf$geneId,
  # pointSize = unit(4, "mm"),
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb"))
)


pkExp_htlist <- anGl_pkExp$an +
  multiProfiles_pkExp$heatmapList + grpHt_pkExp


## row order as decreasing polII signal
if( all(rownames(pkExp_htlist@ht_list[[ polIIData$profileName[1] ]]@matrix) == peakExpDf$gene) ){
  ## polII signal order
  pkExpRowOrd <- order(peakExpDf[[polIICols$exp[1]]], decreasing = TRUE)
  
  ## hasPeak column order
  pkExpRowOrd <- order(peakExpDf[[ tfCols$hasPeak[1] ]],
                       peakExpDf[[ tfCols$hasPeak[2] ]],
                       peakExpDf[[ polIICols$exp[1] ]], decreasing = TRUE)
  
}


# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_pkExp <- paste(comparisonName, ": macs2 peak target genes and top 10% polII signal genes", collapse = "")

# png(filename = paste(outPrefix_pkExp, "_profiles.png", sep = ""), width=wd, height=3500, res = 270)
pdf(file = paste(outPrefix_pkExp, "_profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(pkExp_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_pkExp,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(7, "mm"),
     row_order = pkExpRowOrd,
     padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()



##################################################################################
## plot top 10% expressed genes

## select the genes which are expressed in polII sample
isExpressedDf <- filter_at(.tbl = expressionData, .vars = polIICols$is_expressed,
                           .vars_predicate = any_vars(. == "TRUE"))

multiProfiles_exp <- multi_profile_plots(
  exptInfo = exptData,
  expressionData = isExpressedDf,
  genesToPlot = isExpressedDf$geneId,
  matSource = matrixType,
  matBins = matrixDim,
  clusters = NULL,
  # clusterColor = multiProfiles_all$profileHeatmaps[[1]]$clusterColor,
  profileColors = colorList,
  expressionColor = multiProfiles_all$expressionColor,
  plotExpression = showExpressionHeatmap,
  ylimFraction = ylimList
)


## gene length annotation
anGl_exp <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = isExpressedDf$geneId,
  # pointSize = unit(4, "mm"),
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb"))
)


exp_htlist <- anGl_exp$an +
  multiProfiles_exp$heatmapList
# tssProfiles_exp$heatmapList + tesProfiles_exp$heatmapList

## row order as decreasing polII signal
if( all(rownames(exp_htlist@ht_list[[ polIIData$profileName[1] ]]@matrix) == isExpressedDf$gene) ){
  
  expRowOrd <- order(isExpressedDf[[polIICols$exp[1]]], decreasing = TRUE)
}

# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_exp <- paste(comparisonName, ": top 10% polII signal genes", collapse = "")

# png(filename = paste(outPrefix_expressed, "_profiles.png", sep = ""), width=wd, height=3500, res = 270)
pdf(file = paste(outPrefix_expressed, "_profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(exp_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_exp,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(7, "mm"),
     row_order = expRowOrd,
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()


##################################################################################
# plot SM genes

smDf <- expressionData[which(expressionData$is_SM_gene), ]

multiProfiles_sm <- multi_profile_plots(
  exptInfo = exptData,
  expressionData = smDf,
  genesToPlot = smDf$geneId,
  matSource = matrixType,
  matBins = matrixDim,
  clusters = smDf,
  clusterColor = multiProfiles_all$profileHeatmaps[[1]]$clusterColor,
  profileColors = colorList,
  expressionColor = multiProfiles_all$expressionColor,
  plotExpression = showExpressionHeatmap,
  ylimFraction = ylimList
)

## gene length annotation
anGl_sm <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = smDf$geneId,
  # pointSize = unit(4, "mm"),
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb"))
)

sm_htlist <- anGl_sm$an + multiProfiles_sm$heatmapList

wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_peak= paste(comparisonName, ": SM genes", collapse = "")

# png(filename = paste(outPrefix_sm, "_profiles.png", sep = ""), width=wd, height=3500, res = 270)
pdf(file = paste(outPrefix_sm, "_profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(sm_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(c(5, multiProfiles_all$plotGaps, 5), "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()



##################################################################################





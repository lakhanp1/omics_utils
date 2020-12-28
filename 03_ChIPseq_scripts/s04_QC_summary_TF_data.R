suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))

## 1) peak enrichment distribution
## 2) peak p-value distribution
## 3) peak annotation pie chart
## 4) combined matrix of peak enrichment distribution plots
## 5) combined matrix of peak p-value distribution plots

rm(list = ls())

##################################################################################
analysisName <- "TF_ChIP_summary"
outDir <- here::here("analysis", "02_QC_TF")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")

file_genes <- here::here("data", "reference_data", "AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")

matrixType <- "2kb_summit"
up <- 2000
down <- 2000
body <- 0
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)


smGenes <- suppressMessages(
  AnnotationDbi::select(
    x = orgDb,
    keys = keys(x = orgDb, keytype = "SM_GENE"),
    columns = c("SM_GENE", "SM_CLUSTER", "SM_ID"),
    keytype = "SM_GENE")) %>% 
  dplyr::mutate(geneType = "SM")


##################################################################################

tfSampleList <- readr::read_tsv(file = file_tf_macs2,  comment = "#")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType)

tfInfoList <- purrr::transpose(tfInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


allPlotData <- NULL

peakCountsDf <- tibble::tibble(
  sampleId = character(), peaks_total = numeric(), peaks_pval20 = numeric(),
  peaks_fe3 = numeric(), peaks_pval20_fe3 = numeric()
)


pdf(file = paste(outPrefix, ".macs2.pdf", sep = ""), width = 15, height = 10,
    onefile = TRUE, pointsize = 10)

i <- 3

for (i in 1:nrow(tfInfo)) {
  
  cat(i, ":", tfInfo$sampleId[i], "\n")
  
  peakType <- tfInfo$peakType[i]
  
  ## extract the peak counts to decide best TF replicate
  peaksGr <- rtracklayer::import(con = tfInfo$peakFile[i], format = peakType)
  peakCountsDf <- dplyr::bind_rows(
    peakCountsDf,
    tibble(
      sampleId = tfInfo$sampleId[i],
      peaks_total = length(peaksGr),
      peaks_pval20 = length(which(peaksGr$pValue >= 20)),
      peaks_fe3 = length(which(peaksGr$signalValue > 3)),
      peaks_pval20_fe3 = length(which(peaksGr$pValue >= 20 | peaksGr$signalValue > 3))
    )
  )
  
  backboneGene <- tfInfo$SM_TF[i]

  ## few TFs are mapped to multiple SM clusters. So preparing the list of SM tf data
  smTfInfo <- suppressMessages(
    AnnotationDbi::select(
      x = orgDb,
      keys = na.omit(backboneGene),
      columns = c("SM_GENE", "SM_CLUSTER", "SM_ID"),
      keytype = "GID")) %>%
    dplyr::filter(!is.na(SM_ID))

  clusterGenes <- suppressMessages(
    AnnotationDbi::select(x = orgDb,
                          keys = smTfInfo$SM_ID,
                          columns = c("GID", "SM_GENE"),
                          keytype = "SM_ID")
  )

  genesToMark <- list(
    SM_cluster = unique(clusterGenes$SM_GENE),
    other_SM_genes = setdiff(smGenes$SM_GENE, clusterGenes$SM_GENE)
  )

  chipSummary <- chip_summary(
    sampleId = tfInfo$sampleId[i],
    peakAnnotation = tfInfo$peakAnno[i],
    peakFile = tfInfo$peakFile[i],
    peakFormat = peakType,
    markTargets = genesToMark,
    pointColor = structure(c("red", "blue"), names = names(genesToMark)),
    pointAlpha = structure(c(1, 0.5), names = names(genesToMark))
  )


  plot(chipSummary$figure)

  allPlotData <- dplyr::bind_rows(allPlotData, chipSummary$data)
  
}

dev.off()

bestReplicate <- dplyr::left_join(x = tfInfo, y = peakCountsDf, by = "sampleId") %>% 
  dplyr::select(!!!colnames(peakCountsDf), SM_TF, copyNumber, condition, timePoint, rep) %>% 
  dplyr::group_by(condition, timePoint) %>% 
  dplyr::arrange(desc(peaks_pval20), .by_group = TRUE) %>% 
  dplyr::mutate(
    totalReps = n(),
    bestRep = row_number()
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(sampleId, starts_with("peaks_"), bestRep, totalReps, everything())

readr::write_tsv(x = bestReplicate, path = paste(outPrefix, ".best_replicates.tab", sep = ""))

##################################################################################
## combined summary plot matrix

theme_plot_matrix <- theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    title = element_text(hjust = 0.5, size = 20),
    strip.text = element_text(size = 8, hjust = 0),
    strip.background = element_rect(fill = "white")
  )

facetCols <- 14

## combined summary plot matrix: peak enrichment
gg_all_enrichment <- ggplot(
  data = allPlotData,
  mapping = aes(x = sampleId, y = peakEnrichment)) +
  geom_quasirandom(color = "#bf9a2d") +
  geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 1)) +
  geom_hline(yintercept = 3, color = "blue") +
  labs(title = "masc2 fold enrichment distribution") +
  facet_wrap(facets = ~ sampleId, ncol = facetCols, scales = "free") +
  theme_plot_matrix

# pdf(file = paste(outPrefix, ".macs2_enrichment.pdf", sep = ""), width = 16, height = 16, onefile = TRUE)
png(filename = paste(outPrefix, ".macs2_enrichment.png", sep = ""), width = 15000, height = 9000, res = 300)
print(gg_all_enrichment)
dev.off()

## combined summary plot matrix: peak p-value
gg_all_pval <- ggplot(
  data = allPlotData,
  mapping = aes(x = sampleId, y = peakPval)) +
  geom_quasirandom(color = "#567a0f") +
  geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 1)) +
  geom_hline(yintercept = 20, color = "red") +
  labs(title = "masc2 p-value distribution") +
  facet_wrap(facets = ~ sampleId, ncol = facetCols, scales = "free") +
  theme_plot_matrix

# pdf(file = paste(outPrefix, ".macs2_enrichment.pdf", sep = ""), width = 16, height = 16, onefile = TRUE)
png(filename = paste(outPrefix, ".macs2_pval.png", sep = ""), width = 15000, height = 9000, res = 300)
print(gg_all_pval)
dev.off()


## combined summary plot matrix: peak width
gg_all_width <- ggplot(
  data = allPlotData,
  mapping = aes(x = sampleId, y = peakWidth)) +
  geom_quasirandom(color = "#8472c0") +
  geom_boxplot(width=0.1, fill = NA, outlier.colour = NA, color = alpha("black", 1)) +
  geom_hline(yintercept = 300, color = "red") +
  labs(title = "macs2 peak width distribution") +
  facet_wrap(facets = ~ sampleId, ncol = facetCols, scales = "free") +
  theme_plot_matrix

# pdf(file = paste(outPrefix, ".peak_width.png", sep = ""), width = 16, height = 16, onefile = TRUE)
png(filename = paste(outPrefix, ".peak_width.png", sep = ""), width = 15000, height = 9000, res = 300)
print(gg_all_width)
dev.off()


##################################################################################




##################################################################################
# ## profile plot
# 
# outDir <- dirname(tfInfo$matFile[i])
# outPrefix <- paste(outDir, "/", tfInfo$sampleId[i], ".raw_peaks_summary", sep = "")
# 
# ## control samples to plot alongside TF ChIP sample
# tfControls <- c("AN10300_sCopy_OE_16h_input_FLAG_ChIPMix55_1",
#                 "AN10300_sCopy_OE_16h_input_FLAG_ChIPMix55_2",
#                 "AN0148_sCopy_OE_16h_HA_ChIPMix46_1",
#                 "AN2025_sCopy_OE_16h_HA_ChIPMix36_1")
# 
# 
# ## create profile matrix of 2kb region around peak summit for control samples
# peaksGr <- rtracklayer::import(con = tfInfo$peakFile[i], format = peakType)
# 
# if(length(peaksGr) > 0){
#   if(peakType == "broadPeak"){
#     mcols(peaksGr)$peak <- round(width(peaksGr) / 2)
#   }
#   
#   peakSummitGr <- GenomicRanges::narrow(x = peaksGr,
#                                         start = pmax(peaksGr$peak, 1),
#                                         width = 1)
#   
#   ctrlSampleInfo <- get_sample_information(
#     exptInfoFile = file_exptInfo,
#     samples = tfControls,
#     dataPath = TF_dataPath,
#     profileMatrixSuffix = matrixType)
#   
#   for (ctrlIdx in 1:nrow(ctrlSampleInfo)) {
#     profileMat <- bigwig_profile_matrix(bwFile = ctrlSampleInfo$bwFile[ctrlIdx],
#                                         regions = peakSummitGr,
#                                         signalName = ctrlSampleInfo$profileName[ctrlIdx],
#                                         extend = c(up, down),
#                                         targetName = "summit",
#                                         storeLocal = T,
#                                         localPath = ctrlSampleInfo$matFile[ctrlIdx])
#   }
#   
#   
#   
#   tempSampleInfo <- dplyr::bind_rows(tfInfo[i, ], ctrlSampleInfo) %>% 
#     dplyr::distinct()
#   
#   tfProfileMat <- import_profile_from_file(
#     file = tfInfo$matFile[i],
#     signalName = tfInfo$profileName[i],
#     selectGenes = peakSummitGr$name,
#     up = matrixDim[1], target = matrixDim[2], down = matrixDim[3],
#     targetType = "point", targetName = "summit" 
#   )
#   
#   quantile(tfProfileMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
#   # tfMeanColor <- colorRamp2(quantile(tfProfileMat, c(0.50, 0.995), na.rm = T), c("white", "red"))
#   tfColorList <- sapply(
#     X = tempSampleInfo$sampleId,
#     FUN = function(x){
#       return(colorRamp2(breaks = quantile(tfProfileMat, c(0.50, 0.995), na.rm = T),
#                         colors = c("white", "red")))
#     }
#   )
#   
#   
#   
#   profilePlots <- multi_profile_plots(exptInfo = tempSampleInfo,
#                                       genesToPlot = peakSummitGr$name,
#                                       profileColors = tfColorList,
#                                       targetType = "point", targetName = "summit",
#                                       clusters = NULL,
#                                       showAnnotation = FALSE,
#                                       matBins = matrixDim,
#                                       column_title_gp = gpar(fontsize = 12))
#   
#   
#   rowOrd_peaks <- order(peakSummitGr$signalValue, decreasing = TRUE)
#   
#   # pdf(file = paste(outPrefix, "_profiles.pdf", sep = ""), width = 18, height = 13)
#   png(file = paste(outPrefix, ".profiles.png", sep = ""), width = 3500, height = 2500, res = 250)
#   
#   ht <- draw(
#     profilePlots$heatmapList,
#     main_heatmap = tempSampleInfo$profileName[1],
#     row_order = rowOrd_peaks,
#     column_title = paste(tfInfo$sampleId[i], "binding signal around 2kb region of macs2 peak summit"),
#     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
#     row_sub_title_side = "left",
#     heatmap_legend_side = "bottom",
#     gap = unit(7, "mm"),
#     padding = unit(rep(0.5, times = 4), "cm")
#   )
#   
#   dev.off()
#   
# }
# 
# 








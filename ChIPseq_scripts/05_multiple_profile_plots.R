library(chipmine)

# devtools::install_github(repo = "lakhanp1/chipmine", ref = "dev_v.2")

## file_sampleInfo should have these columns:
## sampleId, bwFile, matFile, IP_tag, profileName, sampleName
file_sampleInfo <- ""
file_genes <- ""

##################################################################################

exptData <- readr::read_tsv(file = file_sampleInfo)

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>%
  dplyr::select(-score) %>%
  dplyr::mutate(length = end - start)

exptData$profileName <- exptData$sampleId
exptData$sampleName <- exptData$sampleId

## for the first time, generate profile matrices.
## next time these matrices can be directly imported
for(i in 1:nrow(exptData)){
  mat <- bigwig_profile_matrix(bwFile = exptData$bwFile[i],
                               regions = file_genes,
                               signalName = exptData$sampleId[i],
                               storeLocal = TRUE, localPath = exptData$matFile[i],
                               extend = c(2000, 1000),
                               target_ratio = 0.4,
                               targetName = "gene")
}


matList <- import_profiles(exptInfo = exptData,
                           geneList = geneSet$geneId,
                           source = "normalizedmatrix",
                           targetType = "point", targetName = "gene",
                           up = 200, target = 200, down = 100)

## tf colors
tfMeanProfile <- NULL
if(length(exptData$sampleId) == 1){
  tfMeanProfile <- matList[[exptData$sampleId[1]]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[exptData$sampleId])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))
tfColorList <- sapply(
  X = exptData$sampleId,
  FUN = function(x){
    return(colorRamp2(breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
                      colors = c("white", "red")))
  }
)

# ylimList <- sapply(X = exptData$sampleId, FUN = function(x) c(0, 80), simplify = FALSE)
ylimList <- list(
  CREEHA_CONTROL4 = c(0, 75), CREEHA_CONTROL5 = c(0, 55), WT_CONTROL5 = c(0, 70),
  CREEHA_10MMAA4 = c(0, 115), CREEHA_10MMAA5 = c(0, 115), WT_10MMAA5 = c(0, 70)
)

##################################################################################

profilePlots <- multi_profile_plots(exptInfo = exptData, genesToPlot = geneSet$geneId,
                                    profileColors = tfColorList,
                                    clusters = NULL,
                                    showAnnotation = FALSE,
                                    matBins = c(200, 200, 100, 10),
                                    matSource = "normalizedmatrix",
                                    column_title_gp = gpar(fontsize = 12),
                                    ylimFraction = ylimList)



# pdf(file = paste(outPrefix, "_profiles.pdf", sep = ""), width = 18, height = 13)
png(file = paste(outPrefix, "_profiles.png", sep = ""), width = 3500, height = 2500, res = 250)

ht <- draw(
  profilePlots$heatmapList,
  main_heatmap = exptData$profileName[1],
  # split = dplyr::select(peakDiffAn, pvalGroup, diffBind),
  column_title = "creE peaks diffbind comparison",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "bottom",
  gap = unit(7, "mm"),
  padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()












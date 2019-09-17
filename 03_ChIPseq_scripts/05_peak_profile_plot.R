library(chipmine)

## this script:
## 1) use diffbind regions and generate profile plot around peak

rm(list = ls())

outDir <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

analysisName <- "creE_diffbind"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

compare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

##################################################################################

file_diffbindInfo <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind", "sampleInfo.txt")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_diffbindRes <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind",
                               "creE_diffbind.annotation.filtered.tab")
TF_dataPath <- here::here("data", "TF_data")

sampleInfo <- suppressMessages(readr::read_tsv(file = file_diffbindInfo))

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   dataPath = TF_dataPath,
                                   profileMatrixSuffix = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

grp1 <- compare[1]
grp2 <- compare[2]
grp1Index <- which(exptData$groupId == grp1)
grp2Index <- which(exptData$groupId == grp2)
grp1Samples <- exptData$sampleId[grp1Index]
grp2Samples <- exptData$sampleId[grp2Index]
grp1SpecificOcc = paste(grp1, ":specific", sep = "")
grp2SpecificOcc = paste(grp2, ":specific", sep = "")
grp1EnrichedCategory <- paste(grp1, ":enriched", sep = "")
grp2EnrichedCategory <- paste(grp2, ":enriched", sep = "")

##################################################################################

# diffDba <- DiffBind::dba.load(file = paste(analysisName, ".dba", sep = ""), dir = outDir, pre = "")

diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffbindRes)) %>% 
  dplyr::select(seqnames, start, end, name, Fold, diffBind, peakOccupancy, categoryDiffbind, pvalFilteredN) %>% 
  dplyr::filter(pvalFilteredN != 0) %>% 
  dplyr::distinct(name, .keep_all = TRUE) %>% 
  dplyr::mutate(
    cluster = dplyr::case_when(
      peakOccupancy == "common" & categoryDiffbind != "common" ~ 
        paste(peakOccupancy, categoryDiffbind, sep = "\n"),
      TRUE ~ categoryDiffbind 
    )
  )

diffGr <- GenomicRanges::makeGRangesFromDataFrame(diffbindRes, keep.extra.columns = T)

#######################
## get the average summit position
peakList <- GenomicRanges::GRangesList(
  lapply(X = exptData$peakFile[c(grp1Index, grp2Index)],
         FUN = rtracklayer::import, format = "narrowPeak")
)

names(peakList) <- exptData$sampleId[c(grp1Index, grp2Index)]

# pgr <- peakList$CREEHA_CONTROL4
## find overlap of each peak GR with diffGr.
## if multiple peaks overlap with a diffGr, select strongest
ovPk <- endoapply(
  X = peakList,
  FUN = function(pgr){
    ovlp <- findOverlaps(query = diffGr, subject = pgr)
    opgr <- pgr[ovlp@to]
    mcols(opgr)$diffGrId <- ovlp@from
    opgr <- as.data.frame(opgr) %>% 
      dplyr::group_by(diffGrId) %>% 
      dplyr::arrange(desc(pValue)) %>% 
      dplyr::slice(1L) %>% 
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    return(opgr)
  })

combinedPeaks <- unlist(ovPk)

summitPos <- GenomicRanges::resize(
  x = GenomicRanges::shift(x = combinedPeaks, shift = combinedPeaks$peak),
  width = 1, fix = "start"
)

avgSummit <- as.data.frame(x = summitPos, row.names = NULL) %>% 
  dplyr::group_by(diffGrId) %>% 
  dplyr::summarise(meanSummit = round(mean(start)))

diffGr$avgSummit[avgSummit$diffGrId] <- avgSummit$meanSummit

#######################

# peakCenterGr <- GenomicRanges::resize(x = diffGr, fix = "center", width = 1, use.names = T)
peakCenterGr <- GRanges(seqnames = seqnames(diffGr),
                        ranges = IRanges(start = diffGr$avgSummit, width = 1))

mcols(peakCenterGr) <- mcols(diffGr)

which(table(diffGr$name) > 1)

peakDiffAn <- dplyr::select(diffbindRes, name, diffBind, peakOccupancy, categoryDiffbind, pvalFilteredN, cluster) %>% 
  dplyr::mutate(pvalGroup = if_else(pvalFilteredN == 0, "weak", "strong")) %>%
  dplyr::rename(geneId = name) %>%
  as.data.frame()

peakDiffAn$diffBind <- factor(x = peakDiffAn$diffBind, levels = c("down", "noDiff", "up"))
peakDiffAn$peakOccupancy <- factor(
  x = peakDiffAn$peakOccupancy,
  levels = c(grp1SpecificOcc, "common", grp2SpecificOcc))
peakDiffAn$cluster <- factor(
  x = peakDiffAn$cluster,
  levels = c(grp1SpecificOcc, paste("common", grp1EnrichedCategory, sep = "\n"),
             "common",
             paste("common", grp2EnrichedCategory, sep = "\n"), grp2SpecificOcc))

# ## for the first time, generate profile matrices.
# ## next time these matrices can be directly imported
# for(i in 1:nrow(exptData)){
#   mat <- bigwig_profile_matrix(bwFile = exptData$bwFile[i],
#                                regions = peakCenterGr,
#                                signalName = exptData$sampleId[i],
#                                storeLocal = TRUE, localPath = exptData$matFile[i],
#                                extend = c(1000, 1000), targetName = "peak center")
# }


matList <- import_profiles(exptInfo = exptData, geneList = diffGr$name,
                           source = "normalizedmatrix",
                           targetType = "point", targetName = "peak center",
                           up = 100, target = 0, down = 100)


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
                      colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ",")))
    )
  }
)

# ylimList <- sapply(X = exptData$sampleId, FUN = function(x) c(0, 80), simplify = FALSE)
ylimList <- list(
  CREEHA_CONTROL4 = c(0, 75), CREEHA_CONTROL5 = c(0, 55), WT_CONTROL5 = c(0, 70),
  CREEHA_10MMAA4 = c(0, 115), CREEHA_10MMAA5 = c(0, 115), WT_10MMAA5 = c(0, 70)
)

profilePlots <- multi_profile_plots(
  exptInfo = exptData, genesToPlot = diffGr$name,
  profileColors = tfColorList,
  clusters = dplyr::select(peakDiffAn, geneId, cluster),
  showAnnotation = FALSE,
  clustOrd = levels(peakDiffAn$cluster),
  targetType = "point",
  targetName = "summit",
  matBins = c(100, 0, 100, 10), matSource = "normalizedmatrix",
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList
)



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






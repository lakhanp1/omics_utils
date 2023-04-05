library(chipmine)
library(ComplexHeatmap)
library(org.AFumigatus.Af293.eg.db)

rm(list = ls())

outDir <- here::here("analysis", "integration_analysis", "diffbind_RNAseq_heatmap")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

##################################################################################
orgDb <- org.AFumigatus.Af293.eg.db

file_diffbindTargets <- here::here("analysis", "ChIPseq_analysis",
                                   "peak_targets", "diffbind_allPeak_targets.tab")
file_diffbindRes <- here::here("analysis", "ChIPseq_analysis", "diffBind", "creE_diffbind.annotation.filtered.tab")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
TF_dataPath <- here::here("data", "TF_data")

chipSamples <- c("CREEHA_CONTROL4", "CREEHA_10MMAA4")
diffbindCompare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

## "CEA17_AA_vs_CEA17_C", "5A9_AA_vs_5A9_C", "5A9_C_vs_CEA17_C", "5A9_AA_vs_CEA17_AA"

# degResultIds <- c("CEA17_AA_vs_CEA17_C", "5A9_AA_vs_5A9_C")
# select_only_degs <- FALSE
analysisName <- "5A9_C_vs_CEA17_C"
degResultIds <- c("5A9_C_vs_CEA17_C")
select_only_degs <- TRUE

outPrefix = paste(outDir, "/", analysisName, ".goodPeaks",  sep = "")

file_deseq2 <- purrr::map_dfr(
  .x = degResultIds,
  .f = function(x){
    list(
      diffPair = x,
      file_diff = here::here("analysis", "RNAseq_data", x, paste(x, ".DEG_all.txt", sep = "")),
      file_deseq2 = here::here("analysis", "RNAseq_data", x, paste(x, ".DESeq2.tab", sep = "")),
      file_rld = here::here("analysis", "RNAseq_data", x, paste(x, ".rlogCounts.tab", sep = ""))
    )
  })

FDR_cut <- 0.05
lfc_cut <- 0.585
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

##################################################################################
grp1 <- diffbindCompare[1]
grp2 <- diffbindCompare[2]
grp1Enrich = paste(grp1, ":enriched", sep = "")
grp2Enrich = paste(grp2, ":enriched", sep = "")
grp1Specific = paste(grp1, ":specific", sep = "")
grp2Specific = paste(grp2, ":specific", sep = "")

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = chipSamples,
                                   dataPath = TF_dataPath,
                                   matrixSource = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered"),
  FUN = function(x){ structure(paste(x, ".", chipSamples, sep = ""), names = chipSamples) },
  simplify = F, USE.NAMES = T
)


## or use org.db
geneSym <- AnnotationDbi::select(x = orgDb,
                                 keys = keys(orgDb, keytype = "GID"),
                                 columns = c("DESCRIPTION"),
                                 keytype = "GID") %>% 
  dplyr::rename(geneId = GID)

##################################################################################
## prepare ChIPseq target data
tf1Id <- chipSamples[1]
tf2Id <- chipSamples[2]

targetSet <- suppressMessages(readr::read_tsv(file = file_diffbindTargets, col_names = T)) %>% 
  dplyr::select(-starts_with("summitSeq."))

targetSet <- dplyr::filter(targetSet, pvalFilteredN > 0)

diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffbindRes, col_names = T)) %>% 
  dplyr::select(name, Fold)

targetSet <- dplyr::left_join(targetSet, diffbindRes, by = "name") %>% 
  dplyr::distinct()
  

dplyr::group_by_at(targetSet, .vars = vars(starts_with("categoryDiffbind"))) %>% 
  dplyr::summarise(n = n())

##################################################################################
## prepare RNAseq data
# 
# diffPairs <- c("CEA17_AA_CEA17_C", "X5A9_AA_CEA17_AA", "X5A9_AA_X5A9_C", "X5A9_C_CEA17_C")
# 
# degData <- suppressMessages(readr::read_tsv(file = file_deg)) %>% 
#   dplyr::filter(Contrast %in% diffPairs)
# 
# degTable <- dplyr::select(degData, Contrast, GeneID, logFC) %>% 
#   tidyr::spread(key = Contrast, value = logFC, fill = 0)
# 
# deseqDf <- suppressMessages(readr::read_tsv(file = file_deseq2)) %>% 
#   dplyr::filter(padj < 0.05)
# 
# degTable <- dplyr::left_join(x = degTable, y = dplyr::select(deseqDf, geneId, log2FoldChange),
#                              by = c("GeneID" = "geneId")) %>% 
#   dplyr::rename(AA_vs_C = log2FoldChange)


## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
get_foldchange <- function(degFile, name){
  
  lfcCol <- "shrinkLog2FC"
  fdrCol <- "padj"
  diffCol <- "diff_shrink_l2fc"
  
  degs <- fread(file = degFile, sep = "\t", header = T, stringsAsFactors = F)
  
  newColName <- structure(c(lfcCol, fdrCol, diffCol),
                          names = paste(c("log2FC.", "padj.", "diff."), name, sep = ""))
  
  
  df <- degs %>%
    dplyr::mutate(
      !!lfcCol := if_else(condition = !!as.name(fdrCol) < FDR_cut,
                          true = !!as.name(lfcCol), false = 0)
    ) %>% 
    tidyr::replace_na(purrr::set_names(list(0), nm = c(lfcCol))) %>% 
    dplyr::select(geneId, !!!c(lfcCol, fdrCol, diffCol)) %>%
    dplyr::rename(!!!newColName )
  
  return(df)
}


i <- 1

for(i in 1:nrow(file_deseq2)){
  dt <- get_foldchange(degFile = file_deseq2$file_diff[i], name = file_deseq2$diffPair[i])
  geneSym <- dplyr::left_join(geneSym, dt, by = c("geneId" = "geneId"))
}

degTable <- dplyr::select(geneSym, geneId, starts_with("log2FC."))
colnames(degTable) <- gsub(pattern = "log2FC.", replacement = "", x = colnames(degTable))

##################################################################################

mergedData <- dplyr::left_join(x = targetSet, y = geneSym, by = c("geneId" = "geneId"))

## optionally select only those genes which are DEGs as per RNAseq
if(select_only_degs){
  mergedData <- dplyr::filter_at(.tbl = mergedData, .vars = vars(starts_with("diff.")),
                                 .vars_predicate = any_vars(. != "noDEG"))
}

readr::write_tsv(x = mergedData, path = paste(outPrefix, ".data.tab", sep = ""))

##################################################################################
## plotting

## select best peak for a gene which has multiple peaks
plotData <- mergedData %>% 
  dplyr::mutate(
    minPeakDist = pmin(abs(!!as.name(tfCols$peakDist[tf1Id])),
                       abs(!!as.name(tfCols$peakDist[tf2Id])),
                       na.rm = TRUE
    ),
    rowKey = paste(name, geneId, sep = ".")
  ) %>% 
  dplyr::group_by(geneId) %>% 
  dplyr::arrange(abs(minPeakDist), desc(bestPval)) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup()

lfcCols <- grep(pattern = "^log2FC.", x = colnames(plotData), perl = T, value = T)

hasPeakMat <- as.matrix(plotData[tfCols$hasPeak[c(tf1Id, tf2Id)]])
row.names(hasPeakMat) <- plotData$rowKey

peakDiffMat <- as.matrix(plotData["Fold"])
row.names(peakDiffMat) <- plotData$rowKey

lfcMat <- as.matrix(plotData[, lfcCols])
row.names(lfcMat) <- plotData$rowKey

plotData$categoryDiffbind <- factor(
  x = plotData$categoryDiffbind,
  levels =  c(grp1Specific, grp1Enrich, "common", grp2Enrich, grp2Specific)
)

ht_opt(
  heatmap_column_names_gp = gpar(fontsize = 16),
  heatmap_row_title_gp = gpar(fontsize = 16),
  legend_title_gp = gpar(fontsize = 16, fontface = "bold"),
  legend_labels_gp = gpar(fontsize = 14)
)

ht1 <- Heatmap(
  matrix = hasPeakMat,
  name = "hasPeak",
  col = colorRamp2(breaks = c(TRUE, FALSE), colors = c("green", "black")),
  column_labels = gsub(pattern = "hasPeak.", replacement = "", colnames(hasPeakMat)),
  cluster_rows = F, cluster_columns = F,
  show_row_names = F,
  heatmap_legend_param = list(
    title = "rglT peak", at = c(TRUE, FALSE), color_bar = "discrete", labels = c("Yes", "No"),
    grid_height = unit(1, "cm"), grid_width = unit(10, "mm")
  ),
  width = unit(2, "cm")
)

ht2 <- Heatmap(
  matrix = peakDiffMat,
  name = "diffbind",
  col = colorRamp2(breaks = c(-5, 0, 5), colors = c("red", "white", "blue")),
  na_col = "white",
  column_labels = c("rglT peak difference\n 10MMAA vs CONTROL"),
  show_row_dend = F, show_column_dend = F, cluster_columns = F,
  show_row_names = F,
  heatmap_legend_param = list(
    title = "\nDiffbind\nlog2(fold-change)",
    at = c(-6, -3, 0, 3, 6),
    legend_height = unit(4, "cm"), grid_width = unit(7, "mm")
  ),
  width = unit(2, "cm")
)

ht3 <- Heatmap(
  matrix = lfcMat,
  name = "lfc",
  col = colorRamp2(breaks = -3:3, colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr")),
  na_col = "white",
  column_labels = gsub(pattern = "log2FC.", replacement = "", colnames(lfcMat)),
  show_row_dend = F, show_column_dend = F, cluster_columns = F,
  show_row_names = F,
  row_title_rot = 0,
  cluster_row_slices = FALSE,
  heatmap_legend_param = list(
    title = "\nRNAseq\nlog2(fold-change)",
    legend_height = unit(4, "cm"), grid_width = unit(7, "mm")
  ),
  width = unit(6, "cm")
)

htList <- ht1 + ht2 + ht3

plotTitle <- paste("rglT ChIPseq and", analysisName, "comparison RNAseq data")

png(file = paste(outPrefix, ".heatmap.png", sep = ""), width = 3000, height = 3000, res = 250)
draw(object = htList,
     main_heatmap = "lfc",
     split = plotData$categoryDiffbind,
     column_title = plotTitle,
     column_title_gp = gpar(fontsize = 20, fontface = "bold"),
     padding = unit(c(10, 2, 2, 2), "mm"),
     row_sub_title_side = "left"
)

dev.off()

pdf(file = paste(outPrefix, ".heatmap.pdf", sep = ""), width = 12, height = 12)
draw(object = htList,
     main_heatmap = "lfc",
     column_title = plotTitle,
     column_title_gp = gpar(fontsize = 20, fontface = "bold"),
     split = plotData$categoryDiffbind,
     padding = unit(c(5, 2, 2, 2), "mm"),
     row_sub_title_side = "left"
)
dev.off()


ht_opt(RESET = TRUE)


##################################################################################



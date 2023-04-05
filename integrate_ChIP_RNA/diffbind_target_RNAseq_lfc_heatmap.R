library(chipmine)
library(ComplexHeatmap)
library(org.AFumigatus.Af293.eg.db)

rm(list = ls())

outDir <- here::here("analysis", "03_integration_analysis", "diffbind_RNAseq_heatmap")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

##################################################################################
orgDb <- org.AFumigatus.Af293.eg.db

file_diffbindTargets <- here::here("analysis", "02_ChIPseq_analysis",
                                   "01_peak_targets", "diffbind_allPeak_targets.tab")
file_diffbindRes <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind",
                               "creE_diffbind.annotation.filtered.tab")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
TF_dataPath <- here::here("data", "TF_data")
file_markGenes <- paste(outDir, "/", "highlight_genes.txt",  sep = "")

chipSamples <- c("CREEHA_CONTROL4", "CREEHA_10MMAA4")
diffbindCompare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

## "CEA17_AA_CEA17_C", "X5A9_AA_X5A9_C", "X5A9_C_CEA17_C", "X5A9_AA_CEA17_AA"

# select_only_degs <- FALSE
analysisName <- "all_edgeR_DEGs"
degResultIds <- c("CEA17_AA_CEA17_C", "X5A9_AA_X5A9_C", "X5A9_C_CEA17_C", "X5A9_AA_CEA17_AA")
select_only_degs <- FALSE

outPrefix <- paste(outDir, "/", analysisName, ".DiffBind_peaks",  sep = "")

file_rnaseq <- here::here("analysis", "01_RNAseq_data", "RNAseq_edgeR_all.txt")

FDR_cut <- 0.05
lfc_cut <- 1          ## 0.585
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
                                   profileMatrixSuffix = "normalizedmatrix")

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
                                 columns = c("GENE_NAME"),
                                 keytype = "GID") %>% 
  dplyr::rename(geneId = GID)

markGenes <- suppressMessages(readr::read_tsv(file = file_markGenes))

##################################################################################
## prepare ChIPseq target data
tf1Id <- chipSamples[1]
tf2Id <- chipSamples[2]

targetSet <- suppressMessages(readr::read_tsv(file = file_diffbindTargets, col_names = T)) %>% 
  dplyr::select(-starts_with("summitSeq."))

targetSet <- dplyr::filter(targetSet, pvalFilteredN > 0)


dplyr::group_by_at(targetSet, .vars = vars(starts_with("categoryDiffbind"))) %>% 
  dplyr::summarise(n = n())

##################################################################################
## prepare RNAseq data

degData <- suppressMessages(readr::read_tsv(file = file_rnaseq)) %>%
  dplyr::filter(Contrast %in% degResultIds) %>% 
  dplyr::filter(FDR <= FDR_cut, abs(logFC) > lfc_cut)

degData$Contrast <- factor(x = degData$Contrast, levels = degResultIds)

degTable <- data.table::dcast(
  data = as.data.table(degData), formula = GeneID ~ Contrast,
  value.var = c("logFC", "FDR", "DECall"),
  sep = "."
) %>% 
  as_tibble()

degCols <- c(
  "GeneID",
  unlist(lapply(X = levels(degData$Contrast), FUN = grep, x = names(degTable), value = TRUE))
)

degTable <- dplyr::select(degTable, degCols)

# ## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
# get_foldchange <- function(degFile, name){
#   
#   lfcCol <- "shrinkLog2FC"
#   fdrCol <- "padj"
#   diffCol <- "diff_shrink_l2fc"
#   
#   degs <- fread(file = degFile, sep = "\t", header = T, stringsAsFactors = F)
#   
#   newColName <- structure(c(lfcCol, fdrCol, diffCol),
#                           names = paste(c("log2FC.", "padj.", "diff."), name, sep = ""))
#   
#   
#   df <- degs %>%
#     dplyr::mutate(
#       !!lfcCol := if_else(condition = !!as.name(fdrCol) < FDR_cut,
#                           true = !!as.name(lfcCol), false = 0)
#     ) %>% 
#     tidyr::replace_na(purrr::set_names(list(0), nm = c(lfcCol))) %>% 
#     dplyr::select(geneId, !!!c(lfcCol, fdrCol, diffCol)) %>%
#     dplyr::rename(!!!newColName )
#   
#   return(df)
# }
# 
# 
# i <- 1
# 
# for(i in 1:nrow(file_deseq2)){
#   dt <- get_foldchange(degFile = file_deseq2$file_diff[i], name = file_deseq2$diffPair[i])
#   geneSym <- dplyr::left_join(geneSym, dt, by = c("geneId" = "geneId"))
# }
# 
# degTable <- dplyr::select(geneSym, geneId, starts_with("log2FC."))
# colnames(degTable) <- gsub(pattern = "log2FC.", replacement = "", x = colnames(degTable))

##################################################################################

mergedData <- dplyr::full_join(x = targetSet, y = degTable, by = c("geneId" = "GeneID")) %>% 
  dplyr::left_join(y = geneSym, by = c("geneId" = "geneId"))

## optionally select only those genes which are DEGs as per RNAseq
if(select_only_degs){
  mergedData <- dplyr::filter_at(.tbl = mergedData, .vars = vars(starts_with("DECall.")),
                                 .vars_predicate = any_vars(. %in% c("up", "down")))
}

readr::write_tsv(x = mergedData, path = paste(outPrefix, ".data.tab", sep = ""))


##################################################################################
## plotting

## select best peak for a gene which has multiple peaks
plotData <- mergedData %>% 
  dplyr::filter(!is.na(DiffBind_region)) %>% 
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

lfcCols <- grep(pattern = "^logFC.", x = colnames(plotData), perl = T, value = T)

plotData <- tidyr::replace_na(
  data = plotData,
  replace = purrr::map(.x = lfcCols, .f = ~ 0) %>% purrr::set_names(nm = lfcCols)
) %>% 
  dplyr::left_join(y = markGenes, by = "geneId")

hasPeakMat <- as.matrix(plotData[tfCols$hasPeak[c(tf1Id, tf2Id)]])
row.names(hasPeakMat) <- plotData$rowKey

peakDiffMat <- as.matrix(plotData["DiffBind.foldChange"])
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

markAn <- rowAnnotation(
  genes = anno_mark(
    at = which(!is.na(plotData$markLabel)),
    labels = na.exclude(plotData$markLabel),
    labels_gp = gpar(fontface = "italic")
  )
)

ht3 <- Heatmap(
  matrix = lfcMat,
  name = "lfc",
  col = colorRamp2(breaks = -3:3, colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr")),
  na_col = "white",
  column_labels = gsub(pattern = "logFC.", replacement = "", colnames(lfcMat)),
  right_annotation = markAn,
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



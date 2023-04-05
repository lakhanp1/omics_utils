library(chipmine)
library(ComplexHeatmap)
library(org.AFumigatus.Af293.eg.db)

rm(list = ls())

outDir <- here::here("analysis", "integration_analysis", "ChIPseq_RNAseq_heatmap")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

##################################################################################
orgDb <- org.AFumigatus.Af293.eg.db

file_targets <- here::here("analysis", "ChIPseq_analysis", "peak_targets", "peak_targets.curated.filtered.tab")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
TF_dataPath <- here::here("data", "TF_data")
chipSamples <- c("CREEHA_CONTROL4", "CREEHA_CONTROL5", "CREEHA_10MMAA4", "CREEHA_10MMAA5")

## "CEA17_AA_vs_CEA17_C", "5A9_AA_vs_5A9_C", "5A9_C_vs_CEA17_C", "5A9_AA_vs_CEA17_AA"

# degResultIds <- c("CEA17_AA_vs_CEA17_C", "5A9_AA_vs_5A9_C")
# select_only_degs <- FALSE

degResultIds <- c("CEA17_AA_vs_CEA17_C", "5A9_AA_vs_5A9_C")
select_only_degs <- TRUE

outPrefix = paste(outDir, "/peak_AA_vs_C_lfc", sep = "")

file_deseq2 <- purrr::map_dfr(
  .x = degResultIds,
  .f = function(x){
    list(
      diffPair = x,
      file_diff = here::here("analysis", "RNAseq_data", x, paste(x, ".DESeq2.tab", sep = "")),
      file_rld = here::here("analysis", "RNAseq_data", x, paste(x, ".rlogCounts.tab", sep = ""))
    )
  })

FDR_cut <- 0.05
lfc_cut <- 0.585
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

##################################################################################
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
s1 <- "CREEHA_CONTROL4"
s2 <- "CREEHA_10MMAA4"

targetSet <- suppressMessages(readr::read_tsv(file = file_targets, col_names = T))

targetSet$groupId <- dplyr::group_indices(
  .data = targetSet,
  !!! lapply(unname(tfCols$hasPeak[c(s1, s2)]), as.name)
)

dplyr::group_by_at(targetSet, .vars = vars(starts_with("hasPeak."), groupId)) %>% 
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
  
  degs <- fread(file = degFile, sep = "\t", header = T, stringsAsFactors = F)
  
  newColName <- structure(c("log2FoldChange", "padj"),
                          names = paste(c("log2FoldChange.", "padj." ), name, sep = ""))
  
  
  df <- degs %>%
    dplyr::mutate(log2FoldChange = if_else(condition = padj < FDR_cut, true = log2FoldChange, false = 0)) %>% 
    tidyr::replace_na(purrr::set_names(list(0), nm = c("log2FoldChange"))) %>% 
    dplyr::select(geneId, log2FoldChange, padj) %>%
    dplyr::rename(!!!newColName )
  
  return(df)
}


i <- 1

for(i in 1:nrow(file_deseq2)){
  dt <- get_foldchange(degFile = file_deseq2$file_diff[i], name = file_deseq2$diffPair[i])
  geneSym <- dplyr::left_join(geneSym, dt, by = c("geneId" = "geneId"))
}

degTable <- dplyr::select(geneSym, geneId, starts_with("log2FoldChange"))
colnames(degTable) <- gsub(pattern = "log2FoldChange.", replacement = "", x = colnames(degTable))

##################################################################################

mergedData <- dplyr::left_join(x = targetSet, y = degTable, by = c("geneId" = "geneId"))

## optionally select only those genes which are DEGs as per RNAseq
if(select_only_degs){
  mergedData <- dplyr::filter_at(.tbl = mergedData, .vars = degResultIds, .vars_predicate = any_vars(. != 0))
}

hasPeakMat <- as.matrix(mergedData[tfCols$hasPeak[c(s1, s2)]])
row.names(hasPeakMat) <- mergedData$geneId

peakPvalMat <- as.matrix(mergedData[tfCols$peakPval[c(s1, s2)]])
row.names(peakPvalMat) <- mergedData$geneId

lfcMat <- as.matrix(mergedData[, degResultIds])
row.names(lfcMat) <- mergedData$geneId



ht1 <- Heatmap(matrix = hasPeakMat,
               name = "hasPeak",
               col = colorRamp2(breaks = c(TRUE, FALSE), colors = c("green", "black")),
               cluster_rows = F, cluster_columns = F,
               show_row_names = F,
               heatmap_legend_param = list(title = "Peak", at = c(TRUE, FALSE), color_bar = "discrete"),
               width = unit(2, "cm")
)

ht2 <- Heatmap(
  matrix = peakPvalMat,
  name = "peakPval",
  col = colorRamp2(breaks = c(1, 50), colors = c("white", "red")),
  na_col = "grey",
  cluster_rows = F, cluster_columns = F,
  show_row_names = F,
  heatmap_legend_param = list(title = "log10(peak-pval)"),
  width = unit(2, "cm")
)

ht3 <- Heatmap(
  matrix = lfcMat,
  name = "lfc",
  col = colorRamp2(breaks = -3:3, colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr")),
  na_col = "white",
  show_row_dend = F, show_column_dend = F, cluster_columns = F,
  show_row_names = F,
  heatmap_legend_param = list(title = "RNAseq log2(fold-change)"),
  width = unit(8, "cm")
)

htList <- ht1 + ht2 + ht3

png(file = paste(outPrefix, "_heatmap.png", sep = ""), width = 800, height = 1000, res = 100)
draw(object = htList,
     main_heatmap = "lfc",
     column_title = "ChIPseq and RNAseq integration for creE data",
     split = mergedData$groupId,
     row_sub_title_side = "left"
)

dev.off()

pdf(file = paste(outPrefix, "_heatmap.pdf", sep = ""), width = 8, height = 10)
draw(object = htList,
     main_heatmap = "lfc",
     column_title = "ChIPseq and RNAseq integration for creE data",
     split = mergedData$groupId,
     row_sub_title_side = "left"
)
dev.off()








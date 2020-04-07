library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(tibble)
library(org.Mmusculus.GRCm38p6.99.eg.db)

## This script
## 1) read the tabular config file for different genesets and plots
## Config file format: TAB delited file
## <deg> <geneId> <title> <output>
## deg: DEG ID from RNAseq_info
## geneId: a ; separated geneIds to be used for plotting
## title: plot title
## output: a suffix for output file
##
## 2) plot1: log2(fold_change) heatmap of specifc genes of interest
## 3) plot2: rlog transformed normalized gene counts heatmap 
## 4) plot1 + plot2 + annotations
## 


rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/topGO_functions.R")
source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

analysisName <- "geneset_plots"
degResult <- "DKO_vs_WT"

diffDataPath <- here::here("analysis", "02_DESeq2_diff")
outDir <- here::here("analysis", "02_DESeq2_diff", degResult, analysisName)

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, "/", degResult, ".geneset", sep = "")

file_geneset <- paste(outDir, "/heatmap.config.tab", sep = "")

file_sampleInfo <- here::here("data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "RNAseq_info.txt")

orgDb <- org.Mmusculus.GRCm38p6.99.eg.db

cutoff_fdr <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc
####################################################################

geneSets <- suppressMessages(readr::read_tsv(file = file_geneset))
rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == degResult)

exptInfo <- read.table(file = file_sampleInfo, header = T, sep = "\t", stringsAsFactors = F)

sampleIds <- unlist(stringr::str_split(string = rnaseqInfo$samples, pattern = ";"))

exptInfo <- dplyr::filter(exptInfo, sampleId %in% sampleIds)


lfcCol <- "log2FoldChange"

## use org.db
geneInfo <- AnnotationDbi::select(
  x = orgDb,
  keys = keys(x = orgDb, keytype = "GID"),
  columns = c("GENE_NAME", "DESCRIPTION"),
  keytype = "GID") %>% 
  dplyr::rename(geneId = GID)

wrap_80 <- scales::wrap_format(80)
####################################################################
## import data
nromCount <- suppressMessages(readr::read_tsv(file = rnaseqInfo$normCount)) %>% 
  dplyr::select(geneId, sampleIds)
fpkmCount <- suppressMessages(readr::read_tsv(file = rnaseqInfo$fpkm)) %>% 
  dplyr::select(geneId, sampleIds)
rldCount <- suppressMessages(readr::read_tsv(file = rnaseqInfo$rld)) %>% 
  dplyr::select(geneId, sampleIds)


## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
get_foldchange <- function(degFile, name, lfcCol = "log2FoldChange", fdrCol = "padj", fdr_cut){
  
  degs <- fread(file = degFile, sep = "\t", header = T, stringsAsFactors = F)
  
  newColName <- structure(c(lfcCol, fdrCol),
                          names = paste(c("lfc.", "padj." ), name, sep = ""))
  
  df <- degs %>%
    dplyr::mutate(!! lfcCol := if_else(condition = !!sym(fdrCol) < fdr_cut,
                                       true = !! as.name(lfcCol), false = 0)) %>% 
    tidyr::replace_na(purrr::set_names(list(0), nm = c(lfcCol))) %>% 
    dplyr::select(geneId, !!lfcCol, !!fdrCol) %>%
    dplyr::distinct() %>% 
    dplyr::rename(!!!newColName )
  
  return(df)
}


i <- 1

for(i in 1:nrow(rnaseqInfo)){
  dt <- get_foldchange(degFile = rnaseqInfo$deseq2[i], name = rnaseqInfo$comparison[i],
                       lfcCol = lfcCol, fdr_cut = cutoff_fdr)
  geneInfo <- dplyr::left_join(geneInfo, dt, by = c("geneId" = "geneId"))
}


geneSets <- dplyr::mutate(
  geneSets,
  geneId = stringr::str_split(string = geneId, pattern = ";")
) %>% 
  dplyr::mutate(
    geneId = purrr::map(geneId, unique)
  )

# genes <- dplyr::distinct(geneSets, geneId, .keep_all = TRUE) %>% 
#   dplyr::left_join(y = geneInfo, by = c("geneId" = "geneId")) %>%
#   dplyr::left_join(y = fpkmCount, by = c("geneId" = "geneId"))
# 
# 
# fwrite(x = genes, file = paste(outPrefix, "_data.tab", sep = ""), sep = "\t", col.names = T, quote = F)
# 
# ## remove the rows with NA values
# ## can add additional filters to select specific genes
# genes <- dplyr::filter_at(.tbl = genes, 
#                           .vars = vars(starts_with("padj.")),
#                           .vars_predicate = all_vars(!is.na(.))
# )
# 

####################################################################


rownameCol <- "GENE_NAME"
showRowNames <- TRUE
rowNameFontSize <- 8
colNameFontSize <- 14


## generate the plot for each row of configuration file
i <- 1

plotData <- tibble::tibble(geneId = unlist(geneSets$geneId[i])) %>% 
  dplyr::left_join(y = geneInfo, by = c("geneId" = "geneId")) %>%
  dplyr::left_join(y = rldCount, by = c("geneId" = "geneId"))

plotTitle <- geneSets$title[i]
plotOutSuffix <- geneSets$output[i]

## fold change heatmap
foldChangeDf <- dplyr::select(plotData, geneId, starts_with("lfc")) %>%
  tibble::column_to_rownames(var = "geneId")

foldChangeMat <- data.matrix(foldChangeDf)

colnames(foldChangeMat) <- gsub(pattern = "lfc\\.(.*)", replacement = "\\1", colnames(foldChangeMat), perl = TRUE)


fcHeatmap <- Heatmap(
  matrix = foldChangeMat,
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"), space = "LAB"),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_labels = plotData$GENE_NAME,
  row_names_gp = gpar(fontsize = rowNameFontSize),
  column_names_gp = gpar(fontsize = colNameFontSize), 
  # width = unit(4, "cm"),
  heatmap_legend_param = list(title = "\nlog2(fold_change)")
)


####################################################################

# ## significant DEG annotation heatmap
# diffAnDf <-  dplyr::select(plotData, geneId, starts_with("diff")) %>%
#   tibble::column_to_rownames(var = "geneId")
# 
# colnames(diffAnDf) <- gsub(pattern = "diff\\((.*)\\)", replacement = "\\1", colnames(diffAnDf), perl = TRUE)
# 
# 
# degAnHeatmap <- Heatmap(matrix = diffAnDf,
#                         col = c("up" = "#fb9a99", "down" = "#a6cee3", "noDEG" = "grey95"),
#                         cluster_rows = FALSE,
#                         cluster_columns = FALSE,
#                         show_row_names = showRowNames,
#                         row_names_gp = gpar(fontsize = rowNameFontSize),
#                         column_names_gp = gpar(fontsize = colNameFontSize), 
#                         width = unit(2, "cm"),
#                         na_col = "grey95",
#                         heatmap_legend_param = list(title = "\nSignificant\nDEG type")
# )


####################################################################

## rld z-score heatmap
rldGeneCounts <- dplyr::select(plotData, geneId, !!!sampleIds) %>%
  tibble::column_to_rownames(var = "geneId")

rldMat <- data.matrix(rldGeneCounts)

rldZscoreMat <- chipmine::scale_matrix_rows(x = rldMat)



## plot main rld score heatmap
rldZscoreHeatmap <- Heatmap(
  rldZscoreMat,
  col = colorRamp2(breaks = c(min(rldZscoreMat), 0, max(rldZscoreMat)),
                   colors = c("green", "black", "red"), space = "LAB"), 
  show_row_names = TRUE,
  row_labels = plotData$GENE_NAME,
  row_names_gp = gpar(fontsize = rowNameFontSize),
  column_names_gp = gpar(fontsize = colNameFontSize), 
  cluster_columns = FALSE, 
  # width = unit(10, "cm"),
  row_names_max_width = unit(15, "cm"), 
  heatmap_legend_param = list(title = "z-score(rld score)", color_bar = "continuous")
) 




####################################################################

# ## log2(fold_change) heatmap with annotation
# htList1 <- fcHeatmap + degAnHeatmap
# 
# 
# # png(filename = paste(outPrefix, ".fc_heatmap.png", sep = ""), width=4000, height=6000, res = 550)
# 
# # pdf(file = paste(outPrefix, ",fc_heatmap.pdf", sep = ""), width = 10, height = 10, onefile = TRUE)
# 
# draw(object = htList1,
#      column_title = plotTitle,
#      row_title = "Genes",
#      column_title_gp = gpar(fontsize = 14)
# )
# 
# # dev.off()


####################################################################

# ## rld heatmap with annotation
# htList2 <- rldZscoreHeatmap + degAnHeatmap
# 
# 
# # png(filename = paste(outPrefix, ".rld_heatmap.png", sep = ""), width=4000, height=6000, res = 550)
# 
# 
# draw(object = htList2,
#      column_title = plotTitle,
#      row_title = "Genes",
#      column_title_gp = gpar(fontsize = 14)
# )
# 
# # dev.off()


####################################################################

## log2(fold_change) + rld heatmap with annotation
htList3 <- fcHeatmap + rldZscoreHeatmap

plotHt <- nrow(plotData) * 0.1 + 1
plotWd <- length(sampleIds)
# png(filename = paste(outPrefix, ".fc_rld_heatmap.png", sep = ""), width=6000, height=6000, res = 550)

pdf(file = paste(outPrefix, ".heatmap_fc_rld.", plotOutSuffix, ".pdf", sep = ""),
    height = 5, width = 10)

draw(
  object = htList3,
  column_title = wrap_80(plotTitle),
  row_title = "Genes",
  column_title_gp = gpar(fontsize = 14)
)

dev.off()


















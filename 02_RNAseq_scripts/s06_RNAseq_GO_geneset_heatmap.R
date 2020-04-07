library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(org.HSapiens.gencodev30.eg.db)
library(here)
library(GO.db)


## This script
## 1) read the tabular config file for different GO term sets
##
## Config file format: TAB delited file
## <deg> <GO> <title> <output>
## deg: DEG ID from RNAseq_info
## GO: a ; separated GO term IDs
## title: plot title
## output: a suffix for output file
##
## 2) generate the heatmap for genes belonging to the GO terms
##



rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/topGO_functions.R")
source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

analysisName <- "geneset_plots"

file_sampleInfo <- here::here("analysis", "09_RNAseq_CTCF", "sample_info_RNAseq.txt")
file_RNAseq_info <- here::here("data", "RNAseq_info.txt", sep = "")
diffDataPath <- here::here("analysis", "09_RNAseq_CTCF", "02_DESeq2_diff")

file_goSet <- paste(diffDataPath, "/RNAseq_GO_geneset.config.tab", sep = "/")

orgDb <- org.HSapiens.gencodev30.eg.db

col_lfc <- "log2FoldChange"
col_fdr <- "padj"
col_geneId <- "ENSEMBL"

cutoff_fdr <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc


###########################################################################
## prepare data
sampleInfo <- suppressMessages(readr::read_tsv(file = file_sampleInfo))

termSet <- suppressMessages(readr::read_tsv(file = file_goSet))


i <- 2


degResult <- termSet$deg[i]

outDir <- file.path(diffDataPath, degResult, analysisName)
outPrefix <- paste(outDir, "/", degResult, ".GO_geneset", sep = "")

if(!dir.exists(outDir)){
  dir.create(outDir)
}

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == degResult)

sampleIds <- unlist(strsplit(x = rnaseqInfo$samples, split = ";"))

diffData <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[1])) 

normCounts <- suppressMessages(readr::read_tsv(file = rnaseqInfo$normCount[1])) %>% 
  dplyr::select(geneId, sampleIds)

rldCounts <- suppressMessages(readr::read_tsv(file = rnaseqInfo$rld[1])) %>% 
  dplyr::select(geneId, sampleIds)



goTerms <- unlist(strsplit(x = termSet$go[i], split = ";"))


geneset <- AnnotationDbi::select(
  x = orgDb, keys = goTerms, columns = c("GID"), keytype = "GOALL"
) %>% 
  dplyr::rename(
    !!col_geneId := GID,
    termId = GOALL
  )

geneset <- dplyr::left_join(x = geneset, y = diffData, by = col_geneId)

readr::write_tsv(x = geneset, path = paste(outPrefix, ".data.tab", sep = ""))


###########################################################################
## volcano plot
termColor <- structure(
  .Data = RColorBrewer::brewer.pal(n = length(unique(geneset$termId)), name = "Set1"),
  names = unique(geneset$termId)
)

termColor <- termColor[!is.na(names(termColor))]

vp <- volcano_plot(df = diffData, fdr_col = col_fdr, lfc_col = col_lfc,
                   ylimit = 40, xlimit = c(-5, 5),
                   fdr_cut = cutoff_fdr, lfc_cut = cutoff_lfc,
                   colorGenes = split(x = geneset$geneId, f = geneset$termId),
                   geneColors = termColor)


vp$plot

pdf(file = paste(outPrefix, ".volcano.pdf", sep = ""), width = 10, height = 10)
vpt
dev.off()


###########################################################################
## heatmap

geneSubset <- dplyr::left_join(x = geneset, y = rldCounts, by = "geneId") %>% 
  dplyr::mutate(id = row_number())

rownameCol <- "GENE_NAME"
showRowNames <- TRUE
rowNameFontSize <- 14
colNameFontSize <- 14


## fold change heatmap
foldChangeDf <- dplyr::select(geneSubset, id, !!col_lfc) %>%
  tibble::column_to_rownames(var = "id")

foldChangeMat <- data.matrix(foldChangeDf)

fcHeatmap <- Heatmap(
  matrix = foldChangeMat,
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"), space = "LAB"),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = geneSubset$GENE_NAME,
  row_names_gp = gpar(fontsize = rowNameFontSize),
  column_names_gp = gpar(fontsize = colNameFontSize), 
  width = unit(2, "cm"),
  heatmap_legend_param = list(title = "\nlog2(fold_change)")
)


geneCounts <- dplyr::select(geneSubset, id, !!!sampleIds) %>%
  column_to_rownames(var = "id")

countMat <- data.matrix(geneCounts)

countZscoreMat <- chipmine::scale_matrix_rows(x = countMat)


## plot main rld score heatmap
rldZscoreHeatmap <- Heatmap(
  matrix = countZscoreMat,
  name = "countHt",
  col = colorRamp2(breaks = c(min(countZscoreMat), 0, max(countZscoreMat)),
                   colors = c("green", "black", "red"), space = "LAB"), 
  row_names_gp = gpar(fontsize = rowNameFontSize),
  column_names_gp = gpar(fontsize = colNameFontSize), 
  cluster_columns = FALSE,
  row_title_rot = 0,
  width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
  heatmap_legend_param = list(title = "z-score(rlog counts)", color_bar = "continuous")
) 


htAn <- rowAnnotation(
  goTerm = geneSubset$termId,
  col = list(goTerm = termColor),
  show_legend = FALSE
)

## log2(fold_change) + rld heatmap with annotation
htList <- htAn + rldZscoreHeatmap + fcHeatmap 


# png(filename = paste(outPrefix, "_fc_rld_heatmap.png", sep = ""), width=6000, height=6000, res = 550)
pdf(file = paste(outPrefix, ".heatmap.pdf", sep = ""), width = 12, height = 14)

draw(object = htList,
     main_heatmap = "countHt",
     column_title = paste("Heatmap of DEGs belonging to GO terms:", analysisName),
     split = geneSubset$termId,
     show_row_dend = FALSE,
     # auto_adjust = FALSE,
     row_sub_title_side = "left",
     column_title_gp = gpar(fontsize = 14)
)

dev.off()





suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(org.GRCm38p6.Ensembl100.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(GO.db))


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
source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/s01_topGO_functions.R")
source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

diffDataPath <- here::here("analysis", "02_DESeq2_diff")
file_sampleInfo <- here::here("data", "reference_data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "reference_data", "DESeq2_DEG_info.txt")

file_goSet <- paste(diffDataPath, "/RNAseq_GO_geneset.config.tab", sep = "/")

orgDb <- org.GRCm38p6.Ensembl100.eg.db

col_lfc <- "log2FoldChange"
col_fdr <- "padj"
col_geneId <- "GID"
col_geneName <- "GENE_NAME"

cutoff_fdr <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc


###########################################################################
## prepare data
sampleInfo <- suppressMessages(readr::read_tsv(file = file_sampleInfo))

termSet <- suppressMessages(readr::read_tsv(file = file_goSet))

setRow <- 2

degResult <- termSet$deg[setRow]

outDir <- paste(diffDataPath, "/", degResult, "/geneset_plots", sep = "")

plotTitle <- paste(degResult, ":", termSet$title[setRow])
outPrefix <- paste(outDir, "/", degResult, ".GO_geneset.", termSet$output[setRow], sep = "")


if(!dir.exists(outDir)){
  dir.create(outDir)
}

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == degResult)

sampleIds <- unlist(strsplit(x = rnaseqInfo$samples, split = ";"))

diffData <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[1])) 

normCounts <- suppressMessages(readr::read_tsv(file = rnaseqInfo$normCount[1])) %>% 
  dplyr::select(geneId, !!!sampleIds)

rldCounts <- suppressMessages(readr::read_tsv(file = rnaseqInfo$rld[1])) %>% 
  dplyr::select(geneId, !!!sampleIds)


goTerms <- unlist(strsplit(x = termSet$go[setRow], split = ";"))


goGenes <- AnnotationDbi::select(
  x = orgDb, keys = goTerms, columns = c("GID"), keytype = "GOALL"
) %>% 
  dplyr::rename(
    geneId = GID,
    termId = GOALL
  )

geneset <- dplyr::left_join(x = goGenes, y = diffData, by = "geneId") %>% 
  dplyr::filter(
    abs(!!sym(col_lfc)) >= cutoff_lfc & !!sym(col_fdr) <= cutoff_fdr
  ) %>% 
  dplyr::mutate(
    !!sym(col_geneName) := if_else(
      condition = is.na(!!sym(col_geneName)), true = geneId, false = !!sym(col_geneName)
    )
  )

readr::write_tsv(x = geneset, path = paste(outPrefix, ".data.tab", sep = ""))


###########################################################################
## volcano plot

if(length(unique(geneset$termId)) <= 9){
  termColor <- structure(
    .Data = RColorBrewer::brewer.pal(n = length(unique(geneset$termId)), name = "Set1"),
    names = unique(geneset$termId)
  )
} else{
  termColor <- base::structure(
    .Data = rainbow(n = length(unique(geneset$termId))),
    names = unique(geneset$termId)
  )
}


termColor <- termColor[!is.na(names(termColor))]

pt_vol <- volcano_plot(
  df = diffData, fdr_col = col_fdr, lfc_col = col_lfc,
  ylimit = 5, xlimit = c(-5, 5),
  fdr_cut = cutoff_fdr, lfc_cut = cutoff_lfc,
  title = plotTitle,
  highlightGenesets = split(x = geneset$geneId, f = geneset$termId),
  genesetColor = termColor
)


pt_vol$plot <- pt_vol$plot +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  ) +
  guides(color=guide_legend(ncol=4))

pdf(file = paste(outPrefix, ".volcano.pdf", sep = ""), width = 10, height = 12)
pt_vol$plot
dev.off()


###########################################################################
## heatmap

geneSubset <- dplyr::left_join(x = geneset, y = rldCounts, by = "geneId") %>% 
  dplyr::mutate(id = row_number())

rownameCol <- col_geneName
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
     column_title = plotTitle,
     split = geneSubset$termId,
     show_row_dend = FALSE,
     # auto_adjust = FALSE,
     row_sub_title_side = "left",
     column_title_gp = gpar(fontsize = 14)
)

dev.off()



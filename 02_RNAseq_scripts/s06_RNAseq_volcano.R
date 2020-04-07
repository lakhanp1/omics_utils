library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(tidyverse)
library(ggrepel)
library(org.Mmusculus.GRCm38p6.99.eg.db)

## This script
## 1) read the tabular config file for different genesets
##
## Config file format: TAB delited file
## <deg> <geneId> <title> <output>
## deg: DEG ID from RNAseq_info
## geneId: a ; separated geneIds to be highlighted in volcano plot
## title: plot title
## output: a suffix for output file
##
## 2) generate the volcano plot and optionally highlight genes



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

file_geneset <- paste(outDir, "/volcano.config.tab", sep = "")

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


lfcCol <- "log2FoldChange"
fdrCol <- "padj"

## use org.db
geneInfo <- AnnotationDbi::select(
  x = orgDb,
  keys = keys(x = orgDb, keytype = "GID"),
  columns = c("GENE_NAME", "DESCRIPTION"),
  keytype = "GID") %>% 
  dplyr::rename(geneId = GID)

wrap_80 <- scales::wrap_format(80)
####################################################################

diffData <- suppressMessages(
  readr::read_tsv(file = rnaseqInfo$deseq2)
)

diffData <- dplyr::left_join(x = diffData, y = geneInfo, by = "geneId")


geneSets <- dplyr::mutate(
  geneSets,
  geneId = stringr::str_split(string = geneId, pattern = ";")
) %>% 
  dplyr::mutate(
    geneId = purrr::map(geneId, unique)
  ) %>% 
  tidyr::unnest(cols = c(geneId)) %>% 
  dplyr::left_join(y = dplyr::select(geneInfo, geneId, GENE_NAME), by = "geneId") %>% 
  tidyr::nest(geneId = geneId, GENE_NAME = GENE_NAME)



i <- 2

plotTitle <- geneSets$title[i]
plotOutSuffix <- geneSets$output[i]


pt_vol <- volcano_plot(
  df = diffData,
  title = plotTitle,
  fdr_col = fdrCol, lfc_col = lfcCol,
  fdr_cut = cutoff_fdr, lfc_cut = cutoff_lfc,
  markGenes = geneSets$GENE_NAME[[i]]$GENE_NAME,
  geneNameCol = "GENE_NAME",
  ylimit = 8, xlimit = c(-7, 7)
)

# pt_vol$plot

ggsave(
  filename = paste(outPrefix, ".volcano.", plotOutSuffix, ".svg", sep = ""),
  plot = pt_vol$plot,
  device = "svg", width = 10, height = 10, units = "in"
)






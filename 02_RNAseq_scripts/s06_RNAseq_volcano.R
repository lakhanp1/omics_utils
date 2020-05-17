suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

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

####################################################################

diffDataPath <- here::here("analysis", "07_polII_diff")
file_geneset <- here::here("analysis", "07_polII_diff", "geneset_volcano.config.tab")

file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")

orgDb <- org.Anidulans.FGSCA4.eg.db
col_geneId <- "GID"
col_geneName <- "GENE_NAME"

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

col_lfc <- "log2FoldChange"
col_fdr <- "padj"

####################################################################

## use org.db
geneInfo <- AnnotationDbi::select(
  x = orgDb,
  keys = keys(x = orgDb, keytype = col_geneId),
  columns = col_geneName,
  keytype = col_geneId) %>% 
  dplyr::rename(geneId = !!sym(col_geneId))


geneSets <- suppressMessages(readr::read_tsv(file = file_geneset))

geneSets <- dplyr::mutate(
  geneSets,
  geneId = stringr::str_split(string = geneId, pattern = ";")
) %>% 
  dplyr::mutate(
    geneId = purrr::map(geneId, unique)
  )

####################################################################

i <- 1

degResult <- geneSets$deg[i]
outDir <- here::here("analysis", "07_polII_diff", degResult, "geneset_plots")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, "/", degResult, ".geneset", sep = "")


rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == geneSets$deg[i])

diffData <- suppressMessages(
  readr::read_tsv(file = rnaseqInfo$deseq2)
)

diffData <- dplyr::left_join(x = diffData, y = geneInfo, by = "geneId")

####################################################################


plotTitle <- paste(degResult, ":", geneSets$title[i])
plotOutSuffix <- geneSets$output[i]


pt_vol <- volcano_plot(
  df = diffData,
  title = plotTitle,
  fdr_col = col_fdr, lfc_col = col_lfc,
  fdr_cut = cutoff_fdr, lfc_cut = cutoff_lfc,
  markGenes = unlist(geneSets$geneId[i]),
  geneNameCol = "GENE_NAME",
  ylimit = 15, xlimit = c(-4, 4)
)


pt_vol$plot <- pt_vol$plot +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  )


png(filename = paste(outPrefix, ".volcano.", plotOutSuffix, ".png", sep = ""), width = 2000, height = 2000, res = 250)
pt_vol$plot
dev.off()

# ggsave(
#   filename = paste(outPrefix, ".volcano.", plotOutSuffix, ".svg", sep = ""),
#   plot = pt_vol$plot,
#   device = "svg", width = 10, height = 10, units = "in"
# )
# 
# 




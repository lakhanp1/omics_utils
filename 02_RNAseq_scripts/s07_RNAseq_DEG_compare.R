library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(here)
library(org.HSapiens.gencodev30.eg.db)


##
## This script plots the 
## 1) compare the two DEG list and plot Upset plot
## 2) perform GO enrichment for each subset of the genes

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")
source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/topGO_functions.R")

analysisName <- "PCL5_vs_Huh7_CTCF_KO"
outDir <- here::here("analysis", "09_RNAseq_CTCF", "03_DEG_compare", analysisName)

if(!dir.exists(outDir)){
  stop("outDir does not exist")
}

outPrefix <- paste(outDir, analysisName, sep = "/")

file_RNAseq_info <- here::here("data", "RNAseq_info.txt")
diffDataPath <- here::here("analysis", "09_RNAseq_CTCF", "02_DESeq2_diff")

degResults <-  c("Huh7_CRI_CTCF_KO_vs_Huh7_CRI_Ctrl",
                 "PLC5_CRI_CTCF_KO_vs_PLC5_CRI_Ctrl")

samples <- c()
plotTitle <- "all DEG comparison"

orgDb <- org.HSapiens.gencodev30.eg.db
file_topGO <- "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/geneid2go.HSapiens.GRCh38p12.topGO.map"

FDR_cut <- 0.05
lfc_cut <- 1
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

####################################################################

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degResults)

lfcCol <- "log2FoldChange"

## use org.db
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(x = orgDb, keytype = "ENSEMBL_VERSION"),
                                  columns = c("ENSEMBL", "GENE_NAME", "DESCRIPTION"),
                                  keytype = "ENSEMBL_VERSION") %>% 
  dplyr::rename(geneId = ENSEMBL_VERSION)


####################################################################
## import data

## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
get_foldchange <- function(degFile, name, lfcCol = "log2FoldChange", ...){
  
  degs <- suppressMessages(readr::read_tsv(file = degFile)) %>% 
    dplyr::mutate(
      diff = dplyr::case_when(
        !!sym(lfcCol) >= up_cut & padj <= FDR_cut ~ "up",
        !!sym(lfcCol) <= down_cut & padj <= FDR_cut ~ "down",
        TRUE ~ "noDEG"
      ),
      contrast = !!name
    ) %>% 
    dplyr::filter(diff != "noDEG") %>% 
    dplyr::select(geneId, !!lfcCol, padj, diff, contrast, ...)
  
  return(degs)
}


i <- 1
degData <- NULL

for(i in 1:nrow(rnaseqInfo)){
  dt <- get_foldchange(degFile = rnaseqInfo$deg[i], name = rnaseqInfo$comparison[i],
                       lfcCol = lfcCol, "ENSEMBL")
  degData <- dplyr::bind_rows(degData, dt)
}


# geneInfo <- dplyr::left_join(geneInfo, dt, by = c("geneId" = "geneId"))

degData <- dplyr::mutate(degData, group = paste(contrast, "_", diff, sep = ""))

dplyr::group_by(degData, diff, contrast) %>% 
  dplyr::summarise(n = n())

####################################################################

geneList <- split(x = degData$geneId, f = degData$group)

title <- paste("DEG overlap in", paste(degResults, collapse = ", "), "\n", collapse = " ")
wrap_80 <- scales::wrap_format(80)
title <- wrap_80(title)

cm <- make_comb_mat(geneList, mode = "distinct")

set_size(cm)
comb_name(cm)
comb_size(cm)
comb_degree(cm)

pt <- UpSet(
  m = cm,
  pt_size = unit(7, "mm"), lwd = 3,
  set_order = names(geneList),
  top_annotation = HeatmapAnnotation(
    foo = anno_empty(border = FALSE),
    "combSize" = anno_text(
      x = paste("(", comb_size(cm), ")", sep = ""),
      just = "center", rot = 0
    ),
    "Intersection\nsize" = anno_barplot(
      x = comb_size(cm), 
      border = FALSE, 
      gp = gpar(fill = "black"), 
      height = unit(2, "cm")
    ),
    annotation_name_side = "left", 
    annotation_name_rot = 0
  ),
  right_annotation = upset_right_annotation(
    m = cm, bar_width = 0.5
  ),
  row_names_max_width = max_text_width(
    set_name(cm), gp = gpar(fontsize = 12)
  ),
  width = unit(12, "cm"), height = unit(6, "cm")
)


pdf(file = paste(outPrefix, ".overlap_upset.pdf", sep = ""), width = 12, height = 8)
# png(filename = paste(outPrefix, ".overlap_upset.png", sep = ""), height = 1500, width = 4000, res = 200)
draw(
  pt,
  column_title = title,
)
dev.off()


####################################################################
## GO enrichment for individual groups
degGroupsDf <- degData %>% 
  pivot_wider(
    id_cols = c(geneId, ENSEMBL),
    names_from = c(contrast),
    values_from = c(diff, log2FoldChange),
    names_sep = "."
  )


degGroupsGo <- dplyr::group_by_at(degGroupsDf, .vars = vars(starts_with("diff."))) %>% 
  dplyr::do(
    topGO_enrichment(goMapFile = file_topGO,
                     genes = .$ENSEMBL,
                     type = "BP",
                     goNodeSize = 5)
  )


readr::write_tsv(x = degGroupsGo, path = paste(outPrefix, ".DEG_overlap.topGO.tab", sep = ""))

####################################################################
## plot GO enrichment figures
degGroupsGo <- suppressMessages(
  readr::read_tsv(file = paste(outPrefix, ".DEG_overlap.topGO.tab", sep = ""))
)


goBarPlots <- degGroupsGo %>% 
  dplyr::group_by_at(.vars = vars(starts_with("diff."))) %>% 
  dplyr::slice(1:10) %>% 
  dplyr::mutate(
    groupId = dplyr::group_indices()
  ) %>% 
  tidyr::unite(col = "group", groupId, starts_with("diff."), remove = FALSE) %>% 
  dplyr::do(
    plots = enrichment_bar(
      df = ., title = paste(title, ":", unique(.$group)),
      pvalCol = "weightedFisher",
      termCol = "Term",
      colorCol = "group",
      countCol = "Significant"
    ),
    group = unique(.$group)
  ) %>% 
  tidyr::unnest(group)


goPlotList <- goBarPlots$plots %>% 
  purrr::set_names(goBarPlots$group)



goPlotList <- cowplot::align_plots(plotlist = goPlotList, align = "v", axis = "rl")

pdf(file = paste(outPrefix, ".DEG_overlap.topGO.barplot.pdf", sep = ""),
    width = 10, height = 4, onefile = TRUE)

for (i in names(goPlotList)) {
  plot(goPlotList[[i]])
}

dev.off()


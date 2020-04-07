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
## 1) create a combined gene list for all DEG comparisons
## 2) perform GO analysis over this combined list
## 3) prepare a table showing richness of each GO term in each DEG list


rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")
source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/topGO_functions.R")

analysisName <- "all_DEGs_compare"
outDir <- here::here("analysis", "09_RNAseq_CTCF", "03_DEG_compare", analysisName)

if(!dir.exists(outDir)){
  stop("outDir does not exist")
}

outPrefix <- paste(outDir, analysisName, sep = "/")

file_RNAseq_info <- here::here("data", "RNAseq_info.txt")
diffDataPath <- here::here("analysis", "09_RNAseq_CTCF", "02_DESeq2_diff")

degResults <-  c("PHH_SB_CTCF_OE_vs_PHH_SB_Ctrl",
                 "PLC5_SB_CTCF_OE_vs_PLC5_SB_Ctrl",
                 "Huh7_CRI_CTCF_KO_vs_Huh7_CRI_Ctrl",
                 "PLC5_CRI_CTCF_KO_vs_PLC5_CRI_Ctrl"
)

samples <- c()
plotTitle <- "all DEG comparison"

orgDb <- org.HSapiens.gencodev30.eg.db
file_topGO <- "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/geneid2go.HSapiens.GRCh38p12.topGO.map"

cutoff_fdr <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

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
        !!sym(lfcCol) >= cutoff_up & padj <= cutoff_fdr ~ "up",
        !!sym(lfcCol) <= cutoff_down & padj <= cutoff_fdr ~ "down",
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

topgoRes <- topGO_enrichment(goMapFile = file_topGO,
                             genes = unique(degData$ENSEMBL),
                             type = "BP",
                             goNodeSize = 5)

readr::write_tsv(x = topgoRes, path = paste(outPrefix, ".combined_DEGs.topGO.tab", sep = ""))

goMapPerGroup <- dplyr::group_by(degData, diff, contrast) %>% 
  dplyr::do(GO_map(
    genes = .$ENSEMBL, goTerms = topgoRes$GO.ID, orgDb = orgDb, keytype = "ENSEMBL"
  )) %>% 
  dplyr::ungroup()


readr::write_tsv(x = goMapPerGroup, path = paste(outPrefix, ".GO_map.tab", sep = ""))

goMapPerGroup <- dplyr::mutate(goMapPerGroup,
                               group = paste(contrast, ":", diff, sep = ""))


goMapTable <- tidyr::pivot_wider(
  data = goMapPerGroup,
  id_cols = c(GOID, TERM, backgroundRatio),
  names_from = c(contrast, diff),
  values_from = c(enrichment),
  names_sep = ":"
)

readr::write_tsv(x = goMapTable, path = paste(outPrefix, ".GO_matrix.tab", sep = ""))

enrichmentMat <- goMapTable %>% 
  dplyr::select(GOID, ends_with("up"), ends_with("down")) %>% 
  tibble::column_to_rownames(var = "GOID") %>% 
  as.matrix()


scalledMat <- chipmine::scale_matrix_rows(x = enrichmentMat, add_attr = FALSE)
nanRows <- which(apply(X = scalledMat, MARGIN = 1, FUN = function(x){any(x %in% c("NaN", "NA", "Inf")) }))

enrichmentMat <- enrichmentMat[-nanRows, ]
scalledMat <- scalledMat[-nanRows, ]


ht <- Heatmap(matrix = scalledMat,
              name = "zscore",
              # col = colorRamp2(breaks = c(0, 0.6), colors = c("white", "red"), space = "LAB"),
              cluster_rows = TRUE,
              clustering_distance_rows = "euclidean",
              cluster_columns = FALSE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 14), 
              # width = unit(2, "cm"),
              heatmap_legend_param = list(title = "\nenrichment")
)

ht2 <- Heatmap(matrix = enrichmentMat,
               name = "enrichment",
               col = colorRamp2(breaks = c(0, 0.6), colors = c("black", "yellow"), space = "LAB"),
               cluster_rows = TRUE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               column_names_gp = gpar(fontsize = 14), 
               # width = unit(2, "cm"),
               heatmap_legend_param = list(title = "\nenrichment")
)


ht + ht2

ht2 + ht











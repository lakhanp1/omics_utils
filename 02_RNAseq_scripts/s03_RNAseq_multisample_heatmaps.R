library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(here)
library(org.HSapiens.gencodev30.eg.db)

##
## This script plots the 
## 1) plot1: log2(fold_change) heatmap of specifc genes of interest
## 2) plot2: rlog transformed normalized gene counts heatmap 
## 3) plot1 + plot2 + annotations
##

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

analysisName <- "PHH_vs_PCL5_CTCF_OE"
outDir <- here::here("analysis", "09_RNAseq_CTCF", "03_DEG_compare", analysisName)

if(!dir.exists(outDir)){
  stop("outDir does not exist")
}

outPrefix <- paste(outDir, analysisName, sep = "/")

file_RNAseq_info <- here::here("data", "RNAseq_info.txt")
diffDataPath <- here::here("analysis", "09_RNAseq_CTCF", "02_DESeq2_diff")

degResults <-  c("PHH_SB_CTCF_OE_vs_PHH_SB_Ctrl",
                 "PLC5_SB_CTCF_OE_vs_PLC5_SB_Ctrl")

samples <- c()
plotTitle <- "all DEG comparison"

orgDb <- org.HSapiens.gencodev30.eg.db

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
get_foldchange <- function(degFile, name, lfcCol = "log2FoldChange"){
  
  degs <- suppressMessages(readr::read_tsv(file = degFile))
  
  newColName <- structure(c(lfcCol, "padj"),
                          names = paste(c("lfc.", "padj." ), name, sep = ""))
  
  df <- degs %>%
    dplyr::mutate(!! lfcCol := if_else(condition = padj < FDR_cut, true = !! as.name(lfcCol), false = 0)) %>% 
    tidyr::replace_na(purrr::set_names(list(0), nm = c(lfcCol))) %>% 
    dplyr::select(geneId, !! lfcCol, padj) %>%
    dplyr::distinct() %>% 
    dplyr::rename(!!!newColName )
  
  return(df)
}


i <- 1

for(i in 1:nrow(rnaseqInfo)){
  dt <- get_foldchange(degFile = rnaseqInfo$deg[i], name = rnaseqInfo$comparison[i],
                       lfcCol = lfcCol)
  geneInfo <- dplyr::left_join(geneInfo, dt, by = c("geneId" = "geneId"))
}


# dplyr::filter_at(.tbl = geneInfo, .vars = vars(starts_with("padj.")), .vars_predicate = all_vars(is.na(.)))
allData <- dplyr::filter_at(.tbl = geneInfo,
                            .vars = vars(starts_with("padj.")),
                            .vars_predicate = any_vars(. < FDR_cut))

####################################################################
rownameCol <- "geneId"
showRowNames <- FALSE
rowNameFontSize <- 14
colNameFontSize <- 14

## fold change heatmap
## fold change heatmap
foldChangeDf <- dplyr::select(allData, !! rownameCol, starts_with("lfc")) %>%
  tibble::column_to_rownames(var = rownameCol)

foldChangeMat <- data.matrix(foldChangeDf)

colnames(foldChangeMat) <- gsub(pattern = "lfc\\.(.*)", replacement = "\\1", colnames(foldChangeMat), perl = TRUE)

diffColor <- colorRamp2(
  breaks = c(-3, -2, -1, -0.75, -0.4, 0, 0.4, 0.75, 1, 2, 3),
  colors = RColorBrewer::brewer.pal(n = 11, name = "PuOr")
)

fcHeatmap <- Heatmap(matrix = foldChangeMat,
                     col = diffColor,
                     cluster_rows = TRUE,
                     clustering_distance_rows = "euclidean",
                     cluster_columns = FALSE,
                     show_row_names = FALSE,
                     row_names_gp = gpar(fontsize = rowNameFontSize),
                     column_names_gp = gpar(fontsize = colNameFontSize), 
                     width = unit(10, "cm"),
                     heatmap_legend_param = list(title = "\nlog2(fold_change)")
)




####################################################################

## log2(fold_change) heatmap with annotation
htList1 <- fcHeatmap

# png(filename = paste(outPrefix, "_fc_heatmap.png", sep = ""), width=4000, height=6000, res = 550)

pdf(file = paste(outPrefix, ".fc_heatmap.pdf", sep = ""), width = 10, height = 10, onefile = TRUE)

draw(object = htList1,
     column_title = plotTitle,
     row_title = "Genes",
     column_title_gp = gpar(fontsize = 14)
)

dev.off()









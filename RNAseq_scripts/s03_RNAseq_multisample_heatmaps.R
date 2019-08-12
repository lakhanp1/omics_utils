library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(here)
library(org.Hsapiens.eg.db)

##
## This script plots the 
## 1) plot1: log2(fold_change) heatmap of specifc genes of interest
## 2) plot2: rlog transformed normalized gene counts heatmap 
## 3) plot1 + plot2 + annotations
##


rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/RNAseq_scripts/02_DESeq2_functions.R")

analysisName <- "combined_deg"
outDir <- here::here("analysis", "03_combined_RNAseq")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, analysisName, sep = "/")

file_sampleInfo <- here::here("data", "sample_info.txt")

degResults <-  c("MHCC97L_A_vs_WT", "MHCC97L_B_vs_WT",
                 "MHCC97L_AB_vs_WT", "MHCC97L_AB_vs_A", "MHCC97L_AB_vs_B")

samples <- c()
plotTitle <- "all DEG comparison"

orgDb <- org.Hsapiens.eg.db

FDR_cut <- 0.05
lfc_cut <- 0
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

####################################################################

diffFiles <- purrr::map_dfr(
  .x = degResults,
  .f = function(x){
    list(
      diffPair = x,
      file_diff = here::here("analysis", "02_DESeq2_diff", x, paste(x, ".DESeq2.tab", sep = "")),
      file_rld = here::here("analysis", "RNAseq_data", x, paste(x, ".rlogCounts.tab", sep = ""))
    )
  })

lfcCol <- "log2FoldChange"
# geneInfo <- readr::read_tsv(file = file_geneInfo) %>%
#   distinct(geneId, .keep_all = T)

## use org.db
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(x = orgDb, keytype = "ENSEMBL_VERSION"),
                                  columns = c("GENE_NAME", "DESCRIPTION"),
                                  keytype = "ENSEMBL_VERSION") %>% 
  dplyr::rename(geneId = ENSEMBL_VERSION)

####################################################################
## import data

## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
get_foldchange <- function(degFile, name){
  
  degs <- fread(file = degFile, sep = "\t", header = T, stringsAsFactors = F)
  
  newColName <- structure(c(lfcCol, "padj"),
                          names = paste(c("log2.", "padj." ), name, sep = ""))
  
  
  df <- degs %>%
    dplyr::mutate(!! lfcCol := if_else(condition = padj < FDR_cut, true = !! as.name(lfcCol), false = 0)) %>% 
    tidyr::replace_na(purrr::set_names(list(0), nm = c(lfcCol))) %>% 
    dplyr::select(geneId, !! lfcCol, padj) %>%
    dplyr::rename(!!!newColName )
  
  return(df)
}


i <- 1

for(i in 1:nrow(diffFiles)){
  dt <- get_foldchange(degFile = diffFiles$file_diff[i], name = diffFiles$diffPair[i])
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
foldChangeDf <- dplyr::select(allData, !! rownameCol, starts_with("log2")) %>%
  tibble::column_to_rownames(var = rownameCol)

foldChangeMat <- data.matrix(foldChangeDf)

colnames(foldChangeMat) <- gsub(pattern = "log2\\.(.*)", replacement = "\\1", colnames(foldChangeMat), perl = TRUE)

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









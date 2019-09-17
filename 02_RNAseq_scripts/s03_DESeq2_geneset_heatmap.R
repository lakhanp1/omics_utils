library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(tibble)
library(org.HSapiens.gencodev30.eg.db)

## This script plots the 
## 1) plot1: log2(fold_change) heatmap of specifc genes of interest
## 2) plot2: rlog transformed normalized gene counts heatmap 
## 3) plot1 + plot2 + annotations
##

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")

analysisName <- "geneset1"
outDir <- here::here("analysis", "04_geneset", analysisName)

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, analysisName, sep = "/")

file_geneList <- paste(outDir, "/genelist.txt", sep = "")

file_sampleInfo <- here::here("data", "sample_info.txt")
file_degInfo <- here::here("analysis", "02_DESeq2_diff", "diff_info.txt")
file_normCounts <- here::here("analysis", "01_count_data", "count_data.normCounts.tab")
file_rld <- here::here("analysis", "01_count_data", "count_data.rlogCounts.tab")
file_fpkm <- here::here("analysis", "01_count_data", "count_data.FPKM.tab")

degResults <-  c("MHCC97L_A_vs_WT", "MHCC97L_B_vs_WT", "MHCC97L_AB_vs_WT",
                 "MHCC97L_AB_vs_A", "MHCC97L_AB_vs_B")

samples <- c("MHCC97L_WT_1", "MHCC97L_WT_2", "MHCC97L_WT_3",
             "MHCC97L_A_1", "MHCC97L_A_2", "MHCC97L_A_3",
             "MHCC97L_B_1", "MHCC97L_B_2", "MHCC97L_B_3",
             "MHCC97L_AB_1", "MHCC97L_AB_2", "MHCC97L_AB_3")

plotTitle <- "geneset 1"

orgDb <- org.HSapiens.gencodev30.eg.db

cutoff_qval <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc
####################################################################

genelist <- suppressMessages(readr::read_tsv(file = file_geneList))
diffInfo <- suppressMessages(readr::read_tsv(file = file_degInfo))
exptInfo <- read.table(file = file_sampleInfo, header = T, sep = "\t", row.names = "sampleId")

lfcCol <- "log2FoldChange"
# geneInfo <- readr::read_tsv(file = file_geneInfo) %>%
#   distinct(geneId, .keep_all = T)

## use org.db
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(x = orgDb, keytype = "ENSEMBL_VERSION"),
                                  columns = c("ENSEMBL", "GENE_NAME", "DESCRIPTION"),
                                  keytype = "ENSEMBL_VERSION") %>% 
  dplyr::rename(geneId = ENSEMBL)


####################################################################
## import data
nromCount <- suppressMessages(readr::read_tsv(file = file_normCounts)) %>% 
  dplyr::select(geneId, !!!samples)
fpkmCount <- suppressMessages(readr::read_tsv(file = file_fpkm)) %>% 
  dplyr::select(geneId, !!!samples)
rldCount <- suppressMessages(readr::read_tsv(file = file_rld)) %>% 
  dplyr::select(geneId, !!!samples)


## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
get_foldchange <- function(degFile, name, lfcCol = "log2FoldChange"){
  
  degs <- fread(file = degFile, sep = "\t", header = T, stringsAsFactors = F)
  
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

for(i in 1:nrow(diffInfo)){
  dt <- get_foldchange(degFile = diffInfo$deseq2[i], name = diffInfo$comparison[i],
                       lfcCol = lfcCol)
  geneInfo <- dplyr::left_join(geneInfo, dt, by = c("ENSEMBL_VERSION" = "geneId"))
}


genes <- dplyr::distinct(genelist, geneId, .keep_all = TRUE) %>% 
  dplyr::left_join(y = geneInfo, by = c("geneId" = "geneId")) %>%
  dplyr::left_join(y = fpkmCount, by = c("ENSEMBL_VERSION" = "geneId"))


fwrite(x = genes, file = paste(outPrefix, "_data.tab", sep = ""), sep = "\t", col.names = T, quote = F)

## remove the rows with NA values
## can add additional filters to select specific genes
genes <- dplyr::filter_at(.tbl = genes, 
                          .vars = vars(starts_with("padj.")),
                          .vars_predicate = all_vars(!is.na(.))
)


####################################################################

rownameCol <- "GENE_NAME"
showRowNames <- TRUE
rowNameFontSize <- 14
colNameFontSize <- 14


## fold change heatmap
foldChangeDf <- dplyr::select(genes, !! rownameCol, starts_with("lfc")) %>%
  tibble::column_to_rownames(var = rownameCol)

foldChangeMat <- data.matrix(foldChangeDf)

colnames(foldChangeMat) <- gsub(pattern = "lfc\\((.*)\\)", replacement = "\\1", colnames(foldChangeMat), perl = TRUE)


fcHeatmap <- Heatmap(matrix = foldChangeMat,
                     col = colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"), space = "LAB"),
                     cluster_rows = TRUE,
                     clustering_distance_rows = "euclidean",
                     cluster_columns = FALSE,
                     show_row_names = FALSE,
                     row_names_gp = gpar(fontsize = rowNameFontSize),
                     column_names_gp = gpar(fontsize = colNameFontSize), 
                     width = unit(4, "cm"),
                     heatmap_legend_param = list(title = "\nlog2(fold_change)")
)


####################################################################

## significant DEG annotation heatmap
diffAnDf <-  dplyr::select(genes, !! rownameCol, starts_with("diff")) %>%
  tibble::column_to_rownames(var = rownameCol)

colnames(diffAnDf) <- gsub(pattern = "diff\\((.*)\\)", replacement = "\\1", colnames(diffAnDf), perl = TRUE)


degAnHeatmap <- Heatmap(matrix = diffAnDf,
                        col = c("up" = "#fb9a99", "down" = "#a6cee3", "noDEG" = "grey95"),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_row_names = showRowNames,
                        row_names_gp = gpar(fontsize = rowNameFontSize),
                        column_names_gp = gpar(fontsize = colNameFontSize), 
                        width = unit(2, "cm"),
                        na_col = "grey95",
                        heatmap_legend_param = list(title = "\nSignificant\nDEG type")
)



####################################################################

## rld z-score heatmap
rldGeneCounts <- dplyr::select(genes, !!rownameCol, !!!samples) %>%
  column_to_rownames(var = rownameCol)

rldMat <- data.matrix(rldGeneCounts)

rldZscoreMat <- chipmine::scale_matrix_rows(x = rldMat)



## plot main rld score heatmap
rldZscoreHeatmap <- Heatmap(
  rldZscoreMat,
  col = colorRamp2(breaks = c(min(rldZscoreMat), 0, max(rldZscoreMat)),
                   colors = c("green", "black", "red"), space = "LAB"), 
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = rowNameFontSize),
  column_names_gp = gpar(fontsize = colNameFontSize), 
  cluster_columns = FALSE, 
  width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
  heatmap_legend_param = list(title = "z-score(rld score)", color_bar = "continuous")
) 




####################################################################

## log2(fold_change) heatmap with annotation
htList1 <- fcHeatmap + degAnHeatmap


# png(filename = paste(outPrefix, "_fc_heatmap.png", sep = ""), width=4000, height=6000, res = 550)

# pdf(file = paste(outPrefix, "_fc_heatmap.pdf", sep = ""), width = 10, height = 10, onefile = TRUE)

draw(object = htList1,
     column_title = plotTitle,
     row_title = "Genes",
     column_title_gp = gpar(fontsize = 14)
)

# dev.off()


####################################################################

## rld heatmap with annotation
htList2 <- rldZscoreHeatmap + degAnHeatmap


# png(filename = paste(outPrefix, "_rld_heatmap.png", sep = ""), width=4000, height=6000, res = 550)


draw(object = htList2,
     column_title = plotTitle,
     row_title = "Genes",
     column_title_gp = gpar(fontsize = 14)
)

# dev.off()


####################################################################

## log2(fold_change) + rld heatmap with annotation
htList3 <- fcHeatmap + rldZscoreHeatmap + degAnHeatmap


png(filename = paste(outPrefix, "_fc_rld_heatmap.png", sep = ""), width=6000, height=6000, res = 550)


draw(object = htList3,
     column_title = plotTitle,
     row_title = "Genes",
     column_title_gp = gpar(fontsize = 14)
)

dev.off()


















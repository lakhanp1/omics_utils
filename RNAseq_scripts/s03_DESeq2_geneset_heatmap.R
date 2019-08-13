library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(tibble)



##
## This script plots the 
## 1) plot1: log2(fold_change) heatmap of specifc genes of interest
## 2) plot2: rlog transformed normalized gene counts heatmap 
## 3) plot1 + plot2 + annotations
##


rm(list = ls())

path <- "G:/Analysis_2/CoreData/14_ZhuBo_RNASeq_ybx1_larvae/Day5_vs_Day5_M_diff/specific_genelist_heatmaps"
setwd(path)

compare <- c("Day5", "Day5_M")

file_geneList <- "genelist_inflamation.txt"

sampleInfoFile <- "G:/Analysis_2/CoreData/14_ZhuBo_RNASeq_ybx1_larvae/Day5_vs_Day5_M_diff/sampleInfo.txt"
geneInfoFile <- "E:/Chris_UM/Database/Zebrafish/GRCz10/ensembl_to_zfin.txt"
rldFile <- "G:/Analysis_2/CoreData/14_ZhuBo_RNASeq_ybx1_larvae/Day5_vs_Day5_M_diff/Day5_vs_Day5_M_rlogCounts.tab"
outPrefix <- "genes_inflamation"


plotTitle <- "inflamation genes"

lfcCol <- "shrinkLog2FC"


##Prepare input data
diffFiles <- data.frame(
  comparison = c("Day5_vs_Day6", "Day5_vs_Day5_M", "Day5_M_vs_Day6"),
  degFiles = c("G:/Analysis_2/CoreData/14_ZhuBo_RNASeq_ybx1_larvae/Day5_vs_Day6_diff/Day5_vs_Day6_DEG_all.txt",
               "G:/Analysis_2/CoreData/14_ZhuBo_RNASeq_ybx1_larvae/Day5_vs_Day5_M_diff/Day5_vs_Day5_M_DEG_all.txt",
               "G:/Analysis_2/CoreData/14_ZhuBo_RNASeq_ybx1_larvae/Day5_M_vs_Day6_diff/Day5_M_vs_Day6_DEG_all.txt"),
  stringsAsFactors = FALSE
)

## remove the comparison in which you are not interested
diffFiles <- dplyr::filter(diffFiles, comparison == "Day5_vs_Day5_M")

geneSym <- fread(file = geneInfoFile, sep = "\t", header = T, stringsAsFactors = F, na.strings = "") %>%
  distinct(geneId, .keep_all = T)


####################################################################

## experiment data
exptInfo <- read.table(file = sampleInfoFile, header = T, sep = "\t", row.names = "sampleId")

## if required, modify the row names
rownames(exptInfo) <- sub("_WT", "_", rownames(exptInfo))

designInfo <- droplevels(subset(exptInfo, condition %in% compare))


## extract the fold change values, DEG status and qval fields for genes of interest
rld <- fread(rldFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)


genes <- fread(file = file_geneList, sep = "\t", header = T, stringsAsFactors = F) %>%
  dplyr::distinct(geneId, .keep_all = TRUE) %>% 
  dplyr::left_join(y = geneSym, by = c("geneId" = "geneId")) %>%
  dplyr::left_join(y = rld, by = c("geneId" = "geneID"))



## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
get_foldchange <- function(df, degFile, name){
  
  degs <- fread(file = degFile, sep = "\t", header = T, stringsAsFactors = F)
  
  newColName <- structure(c(lfcCol, "diff_l2fc", "padj"), names = paste(c("log2(", "diff(", "padj(" ), name, ")", sep = ""))
  
  df <- dplyr::left_join(x = df, y = degs, by = c("geneId" = "geneID")) %>%
    dplyr::select(geneId, !! lfcCol, padj, diff_l2fc) %>%
    dplyr::rename(!!!newColName )
  
  return(df)
}



for(i in 1:nrow(diffFiles)){
  dt <- get_foldchange(df = genes, degFile = diffFiles$degFiles[i], name = diffFiles$comparison[i])
  genes <- dplyr::left_join(genes, dt, by = c("geneId" = "geneId"))
}



fwrite(x = genes, file = paste(outPrefix, "_data.tab", sep = ""), sep = "\t", col.names = T, quote = F)

## remove the rows with NA values
## can add additional filters to select specific genes
genes <- dplyr::filter_at(.tbl = genes, 
                         .vars = vars(starts_with("log2")),
                         .vars_predicate = all_vars(!is.na(.))
)

####################################################################

rownameCol <- "geneName"
showRowNames <- TRUE
rowNameFontSize <- 14
colNameFontSize <- 14


## fold change heatmap
foldChangeDf <- dplyr::select(genes, !! rownameCol, starts_with("log2")) %>%
  tibble::column_to_rownames(var = rownameCol)

foldChangeMat <- data.matrix(foldChangeDf)

colnames(foldChangeMat) <- gsub(pattern = "log2\\((.*)\\)", replacement = "\\1", colnames(foldChangeMat), perl = TRUE)


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
rldGeneCounts <- dplyr::select(genes, !! rownameCol, rownames(designInfo)) %>%
  column_to_rownames(var = rownameCol)

rldMat <- data.matrix(rldGeneCounts)

rldZscoreMat <- matrix(data = NA, nrow = nrow(rldMat), 
                      ncol = ncol(rldMat), 
                      dimnames = list(rownames(rldMat), colnames(rldMat))
)

for(i in 1:nrow(rldMat)){
  rldZscoreMat[i,] <- (rldMat[i,] - mean(rldMat[i,]))/sd(rldMat[i,])
}


## plot main rld score heatmap
rldZscoreHeatmap <- Heatmap(rldZscoreMat,
                            col = colorRamp2(breaks = c(min(rldZscoreMat), 0, max(rldZscoreMat)), c("green", "black", "red"), space = "LAB"), 
                            show_row_names = FALSE,
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


















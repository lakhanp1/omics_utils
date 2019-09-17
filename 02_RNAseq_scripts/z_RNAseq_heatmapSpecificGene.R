library(gplots)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(tibble)


# This script plots the heatmap of specifc genes of interest belonging to different category
rm(list = ls())

path = "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/combinedAnalysis"
setwd(path)

##Prepare input data
compare = c("Day5", "Day5_M", "Day6")

outPrefix = paste(paste(compare, collapse = "_vs_"), "DEGs", sep = "_")

sampleInfoFile = "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/sampleInfo.txt"
rldFile = "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/Day5_Day5_M_diff/Day5_vs_Day5_M_rlogCounts.tab"
normCountFile = "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/Day5_Day5_M_diff/Day5_vs_Day5_M_normCounts.tab"
diffFile = "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/Day5_Day5_M_diff/Day5_vs_Day5_M_DEG_all.txt"

# File with Gene name and category column
specificGenes <- read.table("specificGeneList.tab", header = T, stringsAsFactors = F)


####################################################################

plotTitle = paste("Expression of genes in conditions\n", paste0(compare, collapse = ", "))


exptInfo = read.table(file = sampleInfoFile, header = T, sep = "\t", row.names = "sampleId")

## if required, modify the row names
rownames(exptInfo) = sub("_WT", "_", rownames(exptInfo))


rld = fread(rldFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)
normCounts = fread(normCountFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)
diffData = fread(diffFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)

designInfo = droplevels(subset(exptInfo, condition %in% compare))


rldDiffData = left_join(x = specificGenes, y = rld, by = c("geneId" = "geneID")) %>%
  left_join(y = diffData, by = c("geneId" = "geneID"))

fwrite(x = rldDiffData, file = paste(outPrefix, "_rldDiff.tab", sep = ""), sep = "\t", quote = F, col.names = T)

rldDiffData = rldDiffData %>%
  mutate(id = paste(geneName, geneId, sep = " | ")) %>%
  filter_at(.vars = vars(ends_with("_meanCount")), .vars_predicate = all_vars(. >= 1)) 

anCols = names(specificGenes)[!(names(specificGenes) %in% c("geneId"))]

rldGeneCategory = rldDiffData %>%
  dplyr::select(id, !!!anCols) %>%
  column_to_rownames(var = "id")


#need to convert to matrix to plot the heatmap in case of multiple comparisons
degClasses <- as.matrix(rldGeneCategory)


#######################################################################
## rld heatmap

rldGeneCounts = rldDiffData %>%
  dplyr::select(id, rownames(designInfo)) %>%
  column_to_rownames(var = "id")

rldMat = data.matrix(rldGeneCounts)

rldZscoreMat = matrix(data = NA, nrow = nrow(rldMat), 
                      ncol = ncol(rldMat), 
                      dimnames = list(rownames(rldMat), colnames(rldMat))
)

for(i in 1:nrow(rldMat)){
  rldZscoreMat[i,] <- (rldMat[i,] - mean(rldMat[i,]))/sd(rldMat[i,])
}



## heatmap annotations: if required
# haRow = HeatmapAnnotation(df = rldGeneCategory, name = "Annotation", which = "row",
#                           width = unit(0.5, "cm"), 
#                           col = list(
#                             Day5_M_vs_Day6 = c("1" = "orange", "0" = "grey"),
#                             Day5_vs_Day5_M = c("1" = "orange", "0" = "grey"),
#                             Day5_vs_Day6 = c("1" = "orange", "0" = "grey")
#                           ), 
#                           show_annotation_name = FALSE, 
#                           annotation_name_gp = gpar(fontsize = 11)
# )


## plot main rld score heatmap
rldZscoreHeatmap <- Heatmap(rldZscoreMat, row_title = "Genes", column_title = plotTitle,
                       col = colorRamp2(breaks = c(min(rldZscoreMat), 0, max(rldZscoreMat)), c("green", "black", "red"), space = "LAB"), 
                       show_row_names = FALSE, row_names_side = "left", 
                       row_names_gp = gpar(fontsize = 5),
                       column_names_gp = gpar(fontsize = 10), 
                       show_row_dend = FALSE, cluster_columns = FALSE, 
                       # split = rldGeneCategory,
                       width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
                       heatmap_legend_param = list(title = "z-score(rld score)", color_bar = "continuous")) 



## plot degClass heatmap (in case of multiple comparisons being plotted)
traceHeatmap <- Heatmap(degClasses, 
                        width = unit(2, "cm"), 
                        cluster_columns = FALSE, show_row_names = FALSE, 
                        column_names_gp = gpar(fontsize = 10), 
                        col = structure(c("#1f78b4", "grey"), names = c(1, 0)), 
                        heatmap_legend_param = list(title = "\nDEG under condition", at = c(1, 0), labels = c("Yes", "No"))
)


png(filename = paste(outPrefix, "_rldDiff.png", sep = ""), width=5000, height=8000, res = 600)

rldZscoreHeatmap + traceHeatmap

dev.off()



#######################################################################
## log2FoldChange heatmap
foldChange = rldDiffData %>%
  dplyr::select(id, log2FoldChange) %>%
  column_to_rownames(var = "id")

l2fcMat = data.matrix(foldChange)

colPal = rev(brewer.pal(11, "Spectral")[c(1,2,3,5,6,7,9,10,11)])

htCol = colorRamp2(breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), colors = colPal, space = "LAB")

# c("red", "red", "black", "green", "green")
htCol2 = colorRamp2(breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2), 
                   colors = colorRampPalette(colors = c("red", "black", "green"))(7), space = "LAB")


rldHeatmap = Heatmap(l2fcMat, 
                     col = htCol,
                     row_title = "Genes",
                     row_title_gp = gpar(fontsize = 15),
                     column_title = plotTitle,
                     column_title_gp = gpar(fontsize = 17),
                     show_row_names = TRUE, row_names_side = "left", 
                     row_names_gp = gpar(fontsize = 8),
                     show_column_names = FALSE,
                     column_names_gp = gpar(fontsize = 10),
                     show_row_dend = FALSE, cluster_columns = FALSE,
                     # split = rldGeneCategory,
                     cluster_rows = FALSE,
                     row_order = order(foldChange$log2FoldChange),
                     width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
                     heatmap_legend_param = list(title = "log2(PV / PG)", 
                                                 color_bar = "continuous",
                                                 legend_height = unit(5, "cm"),
                                                 title_gp = gpar(fontsize = 15), 
                                                 labels_gp = gpar(fontsize = 10))
                     )


png(filename = paste(outPrefix, "_foldChange.png", sep = ""), width=4000, height=5000, res = 400)
# tiff(filename = paste(outPrefix, "_foldChange.tif", sep = ""), width=4000, height=6000, res = 400)

rldHeatmap + haRow

dev.off()

#######################################################################
## normalized count heatmap














########### ***************************************
#Order geneIdToName same as expressionMat : Do this
########### ***************************************









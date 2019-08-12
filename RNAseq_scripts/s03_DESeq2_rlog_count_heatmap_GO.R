library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dplyr)
library(RColorBrewer)
library(dynamicTreeCut)
library(data.table)
library(tibble)


## This script plots heatmap for all DEGs for given pair of samples. It clusters the genes based ward.d method. Then inbetween, you have to try cutting clusters into different N subclusters such that each clustere will have sufficient genes. 
## Then this cluster cut information is written to the file
## After this step, find the highly enriched GO term for each cluster and add a third coulmn in the file written last step
## This updated file is again read and a heatmap annotation for each cluster is drawn

rm(list = ls())

path = "E:/Chris_UM/Analysis/CoreData/33_ZhuBo_RNASeq3/PG_WT_50dpf_vs_PG_WT_8mpf"
setwd(path)

##Prepare input data
compare = c("PG_WT_50dpf", "PG_WT_8mpf")

sampleInfoFile = "sampleInfo.txt"
normCountFile = "PG_WT_50dpf_vs_PG_WT_8mpf_normCounts.tab"
rldFile = "PG_WT_50dpf_vs_PG_WT_8mpf_rlogCounts.tab"
diffFile = "PG_WT_50dpf_vs_PG_WT_8mpf_DEG_all.txt"

outFilePrefix = paste0(compare, collapse = "_vs_")
p_cutoff = 0.05
lfc_cut = 1
up_cut = lfc_cut
down_cut = lfc_cut * -1


# plotTitle <- paste("Top ",topN, " regulated genes for each sample pair comparison", sep = " ")
plotTitle = paste("Heatmap of normalized counts for significant DEGs in\n", paste0(compare, collapse = " vs "), "sample comparison")
# plotTitle = "Heatmap of normalized counts for significant DEGs\n between PG and PV stages"

exptInfo = read.table(file = sampleInfoFile, header = T, sep = "\t", row.names = "sampleId")

## if required, modify the row names
# rownames(exptInfo) = sub("_WT", "_", rownames(exptInfo))


rld = fread(rldFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)
normCounts = fread(normCountFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)
diffData = fread(diffFile, sep = "\t", stringsAsFactors = F, header = T, data.table = F)


## select only those sample rows which are part of current comparison
designInfo = droplevels(subset(exptInfo, condition %in% compare))

degGenes = filter(diffData, padj < p_cutoff, log2FoldChange >= up_cut | log2FoldChange <= down_cut)


degCounts = left_join(x = degGenes, y = rld, by = c("geneID")) %>%
  mutate(id = paste(geneName, geneID, sep = "|")) %>%
  dplyr::select(geneID, rownames(designInfo)) %>%
  column_to_rownames(var = "geneID")

expressionMat = data.matrix(degCounts)

# expressionMat = log2(expressionMat + 1)

# calculate dendrogram from original polII values: trueExpr[,c("GZ792_vs_diploid", "GZ892_vs_diploid")]
dend1 = hclust(dist(expressionMat), method = "ward.D")
plot(dend1)
clusterCut1 = cutree(dend1, 4)
clusterData1 = data.frame(GeneID = names(clusterCut1), cluster = clusterCut1)


# z-score with respect to each gene
# i = 1
zScoreMat = matrix(data = NA, nrow = nrow(expressionMat), 
                   ncol = ncol(expressionMat), 
                   dimnames = list(rownames(expressionMat), colnames(expressionMat))
                   )

for(i in 1:nrow(expressionMat)){
  zScoreMat[i,] <- (expressionMat[i,] - mean(expressionMat[i,]))/sd(expressionMat[i,])
}


# calculate dendrogram from original polII values: trueExpr[,c("GZ792_vs_diploid", "GZ892_vs_diploid")]
dend2 = hclust(dist(expressionMat), method = "ward.D")
plot(dend2)
clusterCut2 = cutree(dend2, 4)
clusterData2 = data.frame(GeneID = names(clusterCut2), cluster = clusterCut2)



ht1 = Heatmap(zScoreMat, column_title = plotTitle,
        col = colorRamp2(breaks = c(min(zScoreMat), 0, max(zScoreMat)), c("green", "black", "red"), space = "LAB"), 
        # top_annotation = haUp, top_annotation_height = unit(3, "cm"), 
        show_row_names = FALSE, row_names_side = "left", column_names_gp = gpar(fontsize = 10),        
        width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
        cluster_rows = TRUE,
        show_row_dend = TRUE, row_dend_width = unit(30, "mm"),
        cluster_columns = TRUE, 
        show_column_dend = TRUE,
        # split = clusterData$cluster,
        row_title = "Significant DEGs",
        heatmap_legend_param = list(title = "\nz-score(rlog counts)", 
                                    color_bar = "continuous",
                                    legend_width = unit(7, "cm"),
                                    legend_height = unit(2, "cm"),
                                    title_gp = gpar(fontsize = 14, lineheight = 0.5, srt = 90),
                                    labels_gp = gpar(fontsize = 12, lineheight = 1.5), 
                                    legend_direction = "vertical")
)


outFileName = paste0(outFilePrefix, "_DEG_heatmap.png", collapse = "")
png(filename = outFileName, width=8000, height=10000, res = 600)

draw(ht1, heatmap_legend_side = "right")

dev.off()

# Dendrogram cluster information is written to the file and the cluster numbers are replaced with highly enriched GO term in each cluster
write.table(clusterData, file = "clusterInfo.temp.tab", quote = F, sep = "\t", row.names = F, col.names = T)

######################################################################################
###################################################################################### 
# 
#                                       STOP HERE
#                          Add the GO term column to the file written 
# 
######################################################################################
######################################################################################


## For each cluster, identify the highly enriched GO. For each gene, write appropriate GO term depending on which cluster it falls into

clusterWithGO = read.table(file = "clusterInfo.temp.tab", header = T, sep = "\t")

rownames(clusterWithGO) = clusterWithGO$GeneID

#heatmap annotation boxplot
haUp = HeatmapAnnotation(b1 = anno_boxplot(zScoreMat, which = "column", axis = TRUE), annotation_height = 10)

#heatmap row annotation
rowAnnoColor = structure(brewer.pal(n = 7, name = "Set1"), names = as.character(unique(clusterWithGO$GO)))
haRow = HeatmapAnnotation(GO_term = clusterWithGO$GO, which = "row", width = unit(0.5, "cm"),
                           col = list(GO_term = rowAnnoColor),
                           annotation_legend_param = list(title = "Enriched GO BP term in each cluster",
                                                          title_gp = gpar(fontsize = 14, fontface = "bold"), 
                                                          labels_gp = gpar(fontsize = 14, lineheight = 0.5),
                                                          grid_height = unit(0.7, "cm"),
                                                          grid_width = unit(0.7, "cm"))
                           )


#plot heatmap
rpkmHeatmap = Heatmap(zScoreMat,
                       col = colorRamp2(breaks = c(min(zScoreMat), 0, max(zScoreMat)), c("green", "black", "red"), space = "LAB"), 
                       top_annotation = haUp, top_annotation_height = unit(3, "cm"), 
                       column_title = plotTitle,
                       column_names_gp = gpar(fontsize = 16, fontface = "bold"),
                       show_row_names = FALSE,
                       show_row_dend = TRUE, cluster_columns = FALSE, row_dend_width = unit(30, "mm"),
                       width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
                       heatmap_legend_param = list(title = "\nz-score of log2(normalized count + 1)", 
                                                   color_bar = "continuous",
                                                   legend_width = unit(7, "cm"),
                                                   legend_height = unit(2, "cm"),
                                                   title_gp = gpar(fontsize = 16, fontface = "bold", lineheight = 0.5),
                                                   labels_gp = gpar(fontsize = 16, lineheight = 1.5), 
                                                   legend_direction = "horizontal"), 
                       column_title_gp = gpar(fontsize = 18, fontface = "bold"),
                       split = clusterWithGO$GO, combined_name_fun = NULL
) 


ht_list = haRow + rpkmHeatmap

outFileName = paste0(outFilePrefix, "_DEG_heatmap.png", collapse = "")
png(filename = outFileName, width=8000, height=10000, res = 600)

draw(ht_list, heatmap_legend_side = "right")

dev.off()






library(cummeRbund)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(dynamicTreeCut)

# This script plots heatmap for all DEGs for given pair of samples. It clusters the genes based on perason correlation of FPKM values. Then inbetween, you have to try cutting clusters into different N subclusters such that each clustere will have sufficient genes. 
# Then this cluster cut information is written to the file
# After this step, find the highly enriched GO term for each cluster and add a third coulmn in the file written last step
# This updated file is again read and a heatmap annotation for each cluster is drawn

rm(list = ls())

path <- "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/cuffdiff5"
setwd(path)


#This file has the pairs of sample names. This will be used to draw a second heatmap with trace for which comparison a row is significant
p1 <- c("Day5")
p2 <- c("Day6")
pairs <- data.frame(p1, p2)
pairs$name <- paste(pairs$p1, pairs$p2, sep = " vs ")

#IDs of the samples to plot in the heatmap
sampleIds <- unique(c(p1, p2))

outFileName <- paste(pairs$p1, "vs", pairs$p2, "allGenes.png", sep = "_")


# plotTitle <- paste("Top ",topN, " regulated genes for each sample pair comparison", sep = " ")
plotTitle <- paste("All significantly differentially expressed genes for", pairs$name[1], "sample comparison")


#read all the cufflinks data using cummeRbund package
cuff <- readCufflinks(dir = path)

#sample names
sample.names<-samples(genes(cuff))


#take the top genes
gene.diff<-diffData(object = genes(cuff), x = pairs$p1[1], y = pairs$p2[1])

#filter using p_value and fold_change cutoff
filtered <- subset(x = gene.diff, status == "OK" & q_value < 0.05 & (log2_fold_change > 1 | log2_fold_change < -1))

#get significantly differentially expressed genes
topGenes <- getGenes(object = cuff, geneIdList = filtered$gene_id, sampleIdList = sampleIds)

#FPKM matrix for top N genes
expressionMat <- log2(data.matrix(repFpkmMatrix(topGenes, fullnames = TRUE)) + 0.005)

#rename the row names of expressionMat with just GeneIDs as this will be required to match the geneIds from genesDF 
geneIdToName <- do.call(rbind, strsplit(row.names(expressionMat), split = "|", fixed = TRUE))
colnames(geneIdToName) <- c("GeneName", "GeneID")
rownames(expressionMat) <- geneIdToName[,"GeneID"]


dend = hclust(as.dist(1-cor(t(expressionMat))), method = "complete")

clusterCut <- cutree(dend, 6)
clusterData <- data.frame(GeneID = names(clusterCut), cluster = clusterCut)


# z-score with respect to each gene
for(i in 1:nrow(expressionMat)){
  expressionMat[i,] <- (expressionMat[i,] - mean(expressionMat[i,]))/sd(expressionMat[i,])
}


Heatmap(expressionMat, column_title = plotTitle,
        col = colorRamp2(breaks = c(min(expressionMat), 0, max(expressionMat)), c("green", "black", "red"), space = "LAB"), 
        # top_annotation = haUp, top_annotation_height = unit(3, "cm"), 
        show_row_names = FALSE, row_names_side = "left", column_names_gp = gpar(fontsize = 10), 
        show_row_dend = TRUE, cluster_columns = FALSE, row_dend_width = unit(30, "mm"),
        width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
        heatmap_legend_param = list(title = "z-score of log2(FPKM + 0.005)", color_bar = "continuous", title_gp = gpar(fontsize = 12)),
        # cluster_rows = dend,
        split = clusterData$cluster, combined_name_fun = NULL
        ) 




# Dendrogram cluster information is written to the file and the cluster numbers are replaced with highly enriched GO term in each cluster
write.table(clusterData, file = "clusterInfo.temp.tab", quote = F, sep = "\t", row.names = F, col.names = T)


######################################################################################
# 
# 
# 
#                                       STOP HERE
#                          Add the GO term column to the file written 
# 
# 
# 
######################################################################################

## For each cluster, identify the highly enriched GO. For each gene, write appropriate GO term depending on which cluster it falls into

clusterWithGO <- read.table(file = "clusterInfo.temp.tab", header = T, sep = "\t")

rownames(clusterWithGO) = clusterWithGO$GeneID

#heatmap annotation boxplot
haUp <- HeatmapAnnotation(b1 = anno_boxplot(expressionMat, which = "column", axis = TRUE), annotation_height = 10)

#heatmap row annotation
rowAnnoColor <- structure(brewer.pal(n = 7, name = "Set1"), names = as.character(unique(clusterWithGO$GO)))
haRow <- HeatmapAnnotation(GO_term = clusterWithGO$GO, which = "row", width = unit(0.5, "cm"),
                           col = list(GO_term = rowAnnoColor),
                           annotation_legend_param = list(title = "Enriched GO BP term in each cluster", title_gp = gpar(fontsize = 12)))


#plot heatmap
rpkmHeatmap <- Heatmap(expressionMat, column_title = plotTitle,
                       col = colorRamp2(breaks = c(min(expressionMat), 0, max(expressionMat)), c("green", "black", "red"), space = "LAB"), 
                       top_annotation = haUp, top_annotation_height = unit(3, "cm"), 
                       show_row_names = FALSE, row_names_side = "left", column_names_gp = gpar(fontsize = 10), 
                       show_row_dend = TRUE, cluster_columns = FALSE, row_dend_width = unit(30, "mm"),
                       width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
                       heatmap_legend_param = list(title = "z-score of log2(FPKM + 0.005)", color_bar = "continuous", title_gp = gpar(fontsize = 12)),
                       split = clusterWithGO$GO, combined_name_fun = NULL
                       ) 


ht_list = haRow + rpkmHeatmap

png(filename = outFileName, width=2000, height=3000, res = 180)

draw(ht_list, heatmap_legend_side = "bottom")

dev.off()


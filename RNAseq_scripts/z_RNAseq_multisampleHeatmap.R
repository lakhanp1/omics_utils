library(cummeRbund)
library(gplots)
library(ComplexHeatmap)
library(circlize)


# This script plots the heatmap of log2(FPKM) values of RNA-seq data for multi-sample analysis. It also plots an additional heatmap which shows to which sample comparison each rows belongs.
path <- "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/cuffdiff5"
setwd(path)

####################################################################
#This file has the pairs of sample names. This will be used to draw a second heatmap with trace for which comparison a row is significant
p1 <- c("Day5")
p2 <- c("Day5_M")
pairs <- data.frame(p1, p2)
pairs$name <- paste(pairs$p1, pairs$p2, sep = " vs ")

#IDs of the samples to plot in the heatmap
sampleIds <- unique(c(p1, p2))

#number of top regulated genes to consider
topN <- 250
####################################################################

# plotTitle <- paste("Top ",topN, " regulated genes for each sample pair comparison", sep = " ")
plotTitle <- paste("Top",topN, "regulated genes for", pairs$name[1], "sample comparison")


#read all the cufflinks data using cummeRbund package
cuff <- readCufflinks(dir = path)

#sample names
sample.names<-samples(genes(cuff))



allIDs <- vector(mode = "character") 

genesDF <- data.frame(row.names = seq(1:topN))

#for each comparison, get the top up
for(i in seq(from = 1, to = nrow(pairs), by = 1)){
  pairName <- paste(pairs$p1[i], pairs$p2[i], sep = "_")
  
  #take the top genes
  gene.diff<-diffData(object = genes(cuff), x = pairs$p1[i], y = pairs$p2[i])
  
  #filter using q_value and fold_change cutoff
  filtered <- subset(x = gene.diff, status == "OK" & q_value < 0.05 & (log2_fold_change > 1 | log2_fold_change < -1))
  
  #get the gene IDs of regulated genes in decreasing order of fold_change value
  topRegulated <- filtered$gene_id[order(abs(filtered$log2_fold_change), decreasing = T)][1:topN]
  
  # print(topRegulated)
  
  genesDF[pairName] <- topRegulated
  
  allIDs <- append(x = allIDs, topRegulated)
}

allIDs <- unique(allIDs)

#select top N gene IDs
topGenes <- getGenes(object = cuff, geneIdList = allIDs, sampleIdList = sampleIds)

#FPKM matrix for top N genes
expressionMat <- log2(data.matrix(repFpkmMatrix(topGenes, fullnames = TRUE)) + 0.005)


#rename the row names of expressionMat with just GeneIDs as this will be required to match the geneIds from genesDF 
geneIdToName <- do.call(rbind, strsplit(row.names(expressionMat), split = "|", fixed = TRUE))
colnames(geneIdToName) <- c("GeneName", "GeneID")
rownames(expressionMat) <- geneIdToName[,"GeneID"]

# create the row annotations for each pair of comparison
geneClasses <- data.frame(row.names = row.names(expressionMat))

for(i in names(genesDF)){
  # geneClasses[[i]] <- is.element(allIDs, genesDF[[i]])    #***allIDs will produce wrong output as the order of genes/rows in expressionMat is not same as allIDs
  geneClasses[[i]] <- as.numeric(is.element(row.names(geneClasses), genesDF[[i]]))
  # print(is.element(row.names(geneClasses), genesDF[[i]]))
}


#need to convert to matrix to plot the heatmap
geneClasses <- as.matrix(geneClasses)

# z-score with respect to each gene
for(i in 1:nrow(expressionMat)){
  expressionMat[i,] <- (expressionMat[i,] - mean(expressionMat[i,]))/sd(expressionMat[i,])
}


# use the GeneNames instead of GeneIds in the heatmap
rownames(expressionMat) <- geneIdToName[, "GeneName"]
# rownames(expressionMat) <- geneIdToName[, "GeneName"]


#heatmap annotation boxplot
haUp <- HeatmapAnnotation(b1 = anno_boxplot(expressionMat, which = "column", axis = TRUE), annotation_height = 10)


#plot heatmap
rpkmHeatmap <- Heatmap(expressionMat, row_title = "Genes", column_title = plotTitle,
                       col = colorRamp2(breaks = c(min(expressionMat), 0, max(expressionMat)), c("green", "black", "red"), space = "LAB"), 
                       top_annotation = haUp, top_annotation_height = unit(3, "cm"), 
                       show_row_names = TRUE, row_names_side = "left", column_names_gp = gpar(fontsize = 10), 
                       show_row_dend = FALSE, cluster_columns = FALSE, 
                       width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
                       heatmap_legend_param = list(title = "z-score of log2(FPKM + 1)", color_bar = "continuous")) 

# seq(0,max(expressionMat), length.out = 3)
# c(min(expressionMat), 0, max(expressionMat))
# seq(min(expressionMat),max(expressionMat), length.out = 3)

traceHeatmap <- Heatmap(geneClasses, width = unit(0.7, "cm"), 
                        cluster_columns = FALSE, show_row_names = FALSE, 
                        column_names_gp = gpar(fontsize = 10), 
                        col = structure(c("blue", "gray"), names = c(1, 0)), 
                        heatmap_legend_param = list(title = "Present/absent in top regulated genes", at = c(1, 0), labels = c("Yes", "No")))


# ht_global_opt(RESET = FALSE)


png(filename = "Day5_vs_Day5_M_topDEG.png", width=2000, height=8000, res = 180)

rpkmHeatmap + traceHeatmap

dev.off()








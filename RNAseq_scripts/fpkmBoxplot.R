library(cummeRbund)
library(gplots)
library(ggplot2)
library(circlize)
library(dplyr)
library(tidyr)

# This script plots the histogram of FPKM values for the genes of interest. The gene list should be provided in a separate file. This file can be TAB delimited where each column list down a different category of genes
rm(list = ls())

path <- "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/cuffdiff5"
setwd(path)

genes <- read.table("boxPlotGeneList.txt", header = F)

#read all the cufflinks data using cummeRbund package
cuff <- readCufflinks(dir = path)

#sample names
sample.names<-samples(genes(cuff))

myGenes <- getGenes(cuff, genes$V1)

###########################################################################
# plot title
samplesToPlot <- c("Day5", "Day5_M")

plotTitle <- paste(paste(samplesToPlot, collapse = ", "), ": Intestine development specific genes", sep = ' ')

imageName <- paste(c(samplesToPlot, "intestine_specific.png"), collapse = "_");

###########################################################################

#build a data frame with columns: sampleName, replicate_0, replicate_1, replicate_2, ...
replicates <- replicates(cuff)
reps <- data.frame("sampleName" = replicates$sample_name, "replicateName" = replicates$rep_name, "replicate" = replicates$replicate)
reps <- spread(reps, replicate, replicateName)


#Order the reps as required. This order will decide the sequence of bars in barplot
reps$sampleName <- factor(reps$sampleName, levels = c("Day5", "Day5_M", "Day6"))
reps <- reps[order(reps$sampleName), ]


#Create empty data frame and set the column names.
allFpkms <- data.frame(matrix(nrow = 0, ncol = length(names(reps)[-1]) + 5))
names(allFpkms) <- c(names(reps)[-1], "mean", "ymin", "ymax", "sample", "gene")


for(sp in 1:nrow(reps)){
  
  if(!(reps$sampleName[sp] %in% samplesToPlot)){
    next()
  }
  
  #get the replicate ID for first sample
  repIds <- as.vector(unlist(reps[sp, -1]))
  
  #get the FPKM matrix
  fpkms <- repFpkmMatrix(object = myGenes, fullnames = T, repIdList = repIds)
  
  #change repIds to hold replicate number
  repIds <- names(reps)[-1]
  
  #change the column names: replicateId is replaced with replicateNumber
  names(fpkms) <- repIds
  
  
  fpkms <- log10(fpkms+1)                         #pseudo count to FPKM values
  # fpkms <- fpkms+1
  fpkms$mean <- rowMeans(fpkms[, repIds])         #mean
  fpkms$ymin <- apply(fpkms[, repIds], 1, min)    #ymin
  fpkms$ymax <- apply(fpkms[, repIds], 1, max)    #ymax
  fpkms$sample <- reps$sampleName[sp]             #sample ID
  
  #add gene name columns
  geneIdToName <- as.data.frame(do.call(rbind, strsplit(row.names(fpkms), split = "|", fixed = TRUE)))
  colnames(geneIdToName) <- c("GeneName", "GeneID")
  fpkms$gene <- geneIdToName$GeneName
  
  #Change the row names
  row.names(fpkms) <- paste(row.names(fpkms), fpkms$sample, sep = "|")
  
  # print(head(fpkms))
  allFpkms <- rbind(allFpkms, fpkms)
}

p <- ggplot(allFpkms, aes(x=gene, y=mean, fill=sample)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25) +
  theme_bw() +
  scale_fill_manual(values = c("Day5" = 'wheat', "Day5_M" = '#E69F00', "Day6" = 'lightblue1', "PV_M" ='lightblue4')) +
  labs(title = plotTitle, x = "Genes", y = "log10(FPKM) value") +
  theme(legend.background = element_rect(colour = "black"), plot.title = element_text(size = rel(1.3)), plot.margin=unit(c(1,1,1,1),"cm"), axis.text.x = element_text(angle = 45))

# , axis.text.x = element_text(angle = 45)



png(filename = imageName, width=2000, height=2000, res = 200)
p
dev.off()



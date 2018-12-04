library(cummeRbund)
library(gplots)
library(circlize)
library(dplyr)

# This script plots the volcano plot for of -log10(pvalue) vs log2(foldChange) for given sample pair
path <- "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/cuffdiff5"
setwd(path)

####################################################################

p1 <- c("Day5")
p2 <- c("Day5_M")
pairs <- data.frame(p1, p2)
pairs$name <- paste(pairs$p1, pairs$p2, sep = " vs ")


# plot title
plotTitle <- paste("Volcano plot for", pairs$name[1], "sample comparison")


#read all the cufflinks data using cummeRbund package
cuff <- readCufflinks(dir = path)

#sample names
sample.names<-samples(genes(cuff))


#take the diff data for pair i
gene.diff<-diffData(object = genes(cuff), x = pairs$p1[1], y = pairs$p2[1])

gene.diff$significant <- "NonSignificant"

#q value significant but no large fold change
gene.diff[which(gene.diff$q_value < 0.05 & abs(gene.diff$log2_fold_change) < 1), "significant"] <- "Significant"

#large fold change but not low enough q value
gene.diff[which(gene.diff$q_value > 0.05 & abs(gene.diff$log2_fold_change) > 1), "significant"] <- "FoldChange"

#q value significant and large fold change
gene.diff[which(gene.diff$q_value < 0.05 & abs(gene.diff$log2_fold_change) > 1), "significant"] <- "Significant&FoldChange"

#Drwa Volcano plot
p <- ggplot(data = gene.diff, mapping = aes(x = log2_fold_change, y = -log10(q_value), color = significant)) +
  geom_point(alpha=0.4, size=1.75) + 
  scale_color_manual(values = c("NonSignificant" = "black", "Significant" = "red", "FoldChange" = "orange", "Significant&FoldChange" = "green"), name = "Significance") +
  xlab("log2(fold_change)") + theme_bw() +
  theme(legend.position = c(0.9, 0.7), legend.background = element_rect(colour = "black"), plot.title = element_text(size = rel(1.5)), plot.margin=unit(c(1,1,1,1),"cm")) +
  ggtitle(plotTitle)


imageName = paste(paste("volcano", pairs$p1[1], "vs", pairs$p2[1], sep = "_"), ".png", sep = "");

png(filename = imageName, width=2300, height=2000, res = 200)
p
dev.off()










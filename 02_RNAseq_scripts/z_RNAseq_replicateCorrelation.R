library(cummeRbund)
library(ggplot2)
library(GGally)
# library(grid)
# library(gridExtra)

# This script plots the correlation matrix showing the Pearson correlation between the replicates for each sample. FPKM values are used for the correlation. 

rm(list = ls())

path <- "E:/Chris_UM/Analysis/CoreData/15_ZhuBo_RNASeq2/cuffdiff4"
setwd(path)
# outpath <- dirname(path = path)

pdf(file = "replicateCorrelation.pdf", onefile = TRUE)

replicates <- 3

cuff <- readCufflinks(dir = path)
#sample names
sample.names<-samples(genes(cuff))

gene.rep.matrix<-repFpkmMatrix(genes(cuff), fullnames = T)

gene.count.matrix<-repCountMatrix(genes(cuff), fullnames = T)

# cor(as.matrix(gene.rep.matrix[,1:3]), method = "spearman")
# pairs(as.matrix(log2(gene.rep.matrix[,1:3])))

# gene.rep.df <- as.data.frame(log10(gene.rep.matrix), row.names = row.names(gene.rep.matrix))

gene.rep.df <- as.data.frame(log2(gene.rep.matrix + 1))

pltList <- list()
i <- 1

for(rep in seq(1, ncol(gene.rep.df),replicates)){
  print(c(rep,rep+2))
  end <- rep+replicates-1
  p1 <- ggpairs(gene.rep.df, columns = rep:end,
                lower = list(continuous = "points", combo = "dot"),
                upper = list(continuous = "cor"),
                diag = list(continuous = "densityDiag"),
                axisLabels = 'show',
                title = replicates(cuff)[rep,2])
  print(p1)
  pltList[[i]] <- p1
  i <- i+1
}


dev.off()

write.table(gene.rep.matrix, file = "fpkmMatrix.txt", row.names = T, col.names = T, sep = "\t", quote = F)
write.table(gene.count.matrix, file = "countMatrix.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# p1 <- ggpairs(gene.rep.df, columns = 1:2,
#               lower = list(continuous = "points", combo = "dot"),
#               upper = list(continuous = "cor"),
#               diag = list(continuous = "densityDiag"),
#               axisLabels = 'show',
#               title = replicates(cuff)[1,2])



library("ggplot2")
require("reshape2")
library(grid)
library(tidyr)

rm(list = ls())
path <- "E:/Chris_UM/Analysis/3_cAlbicanceVirulance/analysis2/alleleAssignment/alleleCov_1000/"
setwd(path)

sampleList <- read.table(file = "sample.list", col.names = "sample")
chrFile <- read.table(file = "allele.pairs", sep = "\t")

#read only the name and mean0 columns from the file
mycols <- rep('NULL', 6);
mycols[c(4,6)] <- NA;


chrAB <- data.frame("chr" = character(0), "start" = integer(0), "end" = integer(0), "delMean0" = numeric(0), "sample" = character(0))

# sn <- 2

#iterate over each sample
for(sampleId in sampleList$sample){
  
  #iterate over each chromosome
  for(i in seq(1, nrow(chrFile), by = 1)){
    
    folder <- paste(sampleId, "alleleCov", sep = "_")
    
    # read the data from file
    f1 <- paste(folder, chrFile$V1[i], sep = "/")
    f2 <- paste(folder, chrFile$V2[i], sep = "/")
    
    chrA <- read.table(file = f1, col.names = c("chr", "start", "end", "name", "totalCov", "mean0"), colClasses = mycols)
    chrB <- read.table(file = f2, col.names = c("chr", "start", "end", "name", "totalCov", "mean0"), colClasses = mycols)
    
    
    row.names(chrA) <- chrA$name
    row.names(chrB) <- chrB$name
    
    #split the name column into chr, start, end columns
    chrA <- separate(data = chrA, col = 'name', into = c("chr", "start", "end"), sep = ":|-", convert = TRUE)
    chrB <- separate(data = chrB, col = 'name', into = c("chr", "start", "end"), sep = ":|-", convert = TRUE)
    
    #Append to the chrAB df
    chrAB <- rbind(chrAB, data.frame("chr" = chrA$chr, "start" = chrA$start, "end" = chrA$end, "delMean0" = chrA$mean0-chrB$mean0, "sample" = sampleId, row.names = rownames(chrA)))
    
  }
}

#Function to return the chromosome size in Mb. This is used for axis lables
chrSizeMb <- function(l){
  l <- paste(l/100000, "Mb", sep = "")
  return(l)
}


yLim = max(abs(quantile(chrAB$delMean0, probs = 0.005, names = F)), abs(quantile(chrAB$delMean0, probs = 0.995, names = F)))

chrAB$allele <- ifelse(chrAB$delMean0 > yLim/10, "Allele A", ifelse(chrAB$delMean0 < -yLim/10, "Allele B", "Undetermined"))

p <- ggplot() + 
  geom_rect(data = chrAB, mapping = aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = allele), alpha = 0.2) +
  scale_fill_manual("Allele", values = c("Allele A" = "red", "Allele B" = "green", "Undetermined" = "gray")) + 
  geom_area(data = chrAB, mapping = aes(x = start, y = delMean0), color = "blue", fill = "blue") + 
  scale_x_continuous(labels = chrSizeMb) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 30, face = "bold"), plot.margin = unit(c(1, 1, 1, 1), "cm"), legend.position = "bottom", legend.text = element_text(size=30, face="bold"), legend.title = element_blank()) +
  ggtitle("CA hploid GZ792, GZ892 and WT allele assignment") + ylab("Read depth") +
  coord_cartesian(ylim = c(yLim, -yLim))

# p + facet_wrap(~ chr, ncol = 1)

png(filename = "CA_haploids.png", width=12000, height=600*length(sampleList$sample), res = 400)
p + facet_grid(sample ~ chr, scales = "free_x", space = "free_x")
dev.off()


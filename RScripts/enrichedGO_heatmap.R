require(xlsx)
library(gplots)
library(dplyr)
library(tidyr)
library(XLConnect)
library(ComplexHeatmap)
library(circlize)

# This script craetes the heatmap for enriched GO terms in RNA-seq data. 
# Input: DAVID GO analysis chart in excel with one additional column: Enrichment in Up/Down regulated genes
rm(list = ls())

path <- "E:/Chris_UM/Analysis/CoreData/14_ZhuBo_RNASeq/DAVID_enrichment"
setwd(path)


excel <- loadWorkbook("GO_BP_allCombinations.xlsx")
sheet_names <- getSheets(excel)
names(sheet_names) <- sheet_names

goData <- data.frame("Up_or_Down" = character(0), "Term" = character(0), "PValue" = numeric(0), "Pair" = character(0))

# Only worksheets in this list will be considered to pull the data
includeList <- c("Day5_vs_Day6")

#Read each worksheet and append the data to goData data frame
for(i in sheet_names){
  
  if(! i %in% includeList){
    next()
  }
  
  file <- read.xlsx(file = "GO_BP_allCombinations.xlsx", colIndex = c(1,3,6), sheetName = i)
  file$Pair <- i
  file <- subset(x = file, PValue <= 0.05)
  
  #add rows to data frame goData
  goData = rbind(goData, file)
}

row.names(goData) = c(1:nrow(goData))

goData$PValue <- -log10(goData$PValue)

# merge the Sample and Up/Down column
goData <- unite(data = goData, col = "class", Pair, Up_or_Down, sep = "-")

# spread the table using dplyr package
goData <- spread(data = goData, key = class, value = PValue, fill = 0)

rownames(goData) <- goData$Term

write.table(goData, file = "enriched_GO.txt", sep = "\t", row.names = F, col.names = T, quote = F)

goData <- goData[, !colnames(goData) %in% c("Term")]

goDataMat <- as.matrix(goData)


# heatmap.2(goDataMat, trace = "none", 
#           density.info = 'none', keysize = 2,
#           margins =c(10,25), dendrogram = "row",
#           cexCol = 1,
#           lmat = rbind(c(0,3),c(2,1), c(0,4)), lwid = c(1,4), lhei = c(1,6,1))
# 
# 
# heatmap.2(goDataMat)

#Plot heatmap
ha1 <- Heatmap(goDataMat, 
               col = colorRamp2(breaks = c(0, max(goDataMat, na.rm = T)), c("grey", "red"), space = "LAB"),
               column_title = "Enriched GO terms (p-value <= 0.05) in Day5_vs_Day6 comparison",
               column_title_gp = gpar(fontface = "bold"),
               row_names_max_width = unit(15, "cm"), column_names_max_height = unit(6, "cm"), width = unit(8, "cm"),
               row_dend_width = unit(3, "cm"), show_column_dend = FALSE,
               na_col = "grey", clustering_distance_rows = "euclidean",
               heatmap_legend_param = list(title = "-LOG10(p-value)", color_bar = "continuous", title_position = 'topcenter', legend_direction = 'horizontal', legend_width = unit(5, "cm"), title_gp = gpar(fontsize = 12)))

draw(ha1, heatmap_legend_side = "bottom")

png(filename = "GO_enrichment_Day5_vs_Day6.png", width=3200, height=2000, res = 250)

draw(ha1, heatmap_legend_side = "bottom")

dev.off()












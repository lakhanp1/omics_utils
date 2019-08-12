require(reshape)
library(hashmap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
require(XLConnect)
options(java.parameters = "- Xmx4g")
xlcFreeMemory()

## This script:
## 1. reads the DAVID GO output chart for up and down regulated genes
## 2. plot the significant GO BP terms -log10(p-value) bar plot
## 3. Write the GO chart to excel file

rm(list = ls())

source(file = "E:/Chris_UM/Codes/GO_enrichment/topGO_functions.R")


path = "E:/Chris_UM/Analysis/CoreData/33_ZhuBo_RNASeq3/PG_WT_50dpf_vs_PG_WT_8mpf"
setwd(path)

files = c(up = "PG_WT_50dpf_vs_PG_WT_8mpf_up_DAVID.txt", down = "PG_WT_50dpf_vs_PG_WT_8mpf_down_DAVID.txt")
outPrefix = "PG_WT_50dpf_vs_PG_WT_8mpf"
comparison = "PG_WT_50dpf vs PG_WT_8mpf"
geneInfoFile = "E:/Chris_UM/Database/Zebrafish/GRCz10/ensembl_to_zfin.txt"

filterTerms = c()

excelOut = paste0(outPrefix, "_enrichment.xlsx", collapse = "")
ensembleIdMap = read.table(file = geneInfoFile, sep = "\t", header = T)

idMap = hashmap(keys = as.vector(ensembleIdMap$geneId), values = as.vector(ensembleIdMap$geneName))

# PValue, Benjamini
pvalField = "PValue"
cutoff = 0.05
data = data.frame(up_or_down = character(0), 
                    Term = character(0), 
                    Count = numeric(0), 
                    richFactor = numeric(0), 
                    stringsAsFactors = F)

data[pvalField] = numeric(0)


## This filter will be used to select on GO_BP and significance <= 0.05
# goFilter = list(~Category == 'GOTERM_BP_DIRECT',
#                  interp(as.formula('~col <= val'), col = as.name(pvalField), val = cutoff))

unlink(excelOut, recursive = FALSE, force = FALSE)
exc = loadWorkbook(excelOut , create = TRUE)
xlcFreeMemory()

# i = names(files)[1]
## read GO_chart for Up and Down genes and write the data into single excel file
for(i in names(files)){
  go = read.table(file = files[i], sep = "\t", header = T, stringsAsFactors = F)
  
  wrkSheet = i
  
  #conver the Genes columns to character
  go$Genes = as.character(go$Genes)
  
  #map the ensembl gene id to official gene symbol
  go$Gene_symbol = mapply(function(x){paste(idMap[[unlist(strsplit(as.character(x), split = ", "))]], collapse = ", ")}, go$Genes)
  
  go = dplyr::arrange(go, Category, !!as.name(pvalField)) %>%
    dplyr::mutate(richFactor = Count / Pop.Hits) 
  
  createSheet(exc, name = wrkSheet)
  createFreezePane(exc, sheet = wrkSheet, 2, 2)
  writeWorksheet(object = exc, data = go, sheet = wrkSheet)
  
  
  go$up_or_down = i 
  
  
  # goData = go %>% filter_(.dots = goFilter)
  # Demo:
  # quo(filter(go, Category == "GOTERM_BP_DIRECT", !!(pvalField) <= 0.05) )
  # quo(filter(go, Category == "GOTERM_BP_DIRECT", !!(as.name(pvalField)) <= 0.05) )
  # quo(filter(go, Category == "GOTERM_BP_DIRECT", UQ(pvalField) <= 0.05) )
  # quo(filter(go, Category == "GOTERM_BP_DIRECT", UQ(as.name(pvalField)) <= 0.05) )
  # 
  data = go %>%
    dplyr::select(Category, up_or_down, Term, !!as.name(pvalField), richFactor, Count) %>% 
    bind_rows(data)
  
}




xlcFreeMemory()
saveWorkbook(exc)

log10Field = paste("log10_", pvalField, sep = "")
data = data %>% dplyr::mutate(!!log10Field := -log10(!! as.name(pvalField) ))

if(length(filterTerms) >= 1){
  data = dplyr::filter(.data = data, !(Term %in% filterTerms))
}


goData = data %>% dplyr::filter(Category == "GOTERM_BP_DIRECT", UQ(as.name(pvalField)) <= 0.05) 
goData$Term = gsub("GO:.*~", "", goData$Term, perl = T)

goData = arrange(goData, up_or_down, !!as.name(log10Field))
goData$up_or_down = factor(goData$up_or_down)
goData$Term = factor(goData$Term, levels = unique(goData$Term[order(goData$up_or_down)]))


###########################################################################
## bar chart fo GO BP: Up and Down DEG GOs together
goBar = ggplot(data = goData) +
  geom_bar(mapping = aes_string(x = "Term", y = log10Field, fill = "up_or_down"), stat = "identity") +
  scale_fill_manual(values = c("up" = 'blue', "down" = 'red'),
                    name = "Differential expression\ntype", 
                    breaks = c("up", "down"), labels = c("Up regulated", "Down regulated")) +
  ggtitle(paste("Enriched GO terms (p-value <= 0.05) in", comparison, "comparison")) +
  ylab("-log10(p-value)") + xlab("Biological Process GO Terms") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) + 
  coord_flip()


png(filename = paste0(outPrefix, "_GO_all_bar.png", collapse = ""), width=18000, height=25000, res = 1200)
print(goBar)
dev.off()

###########################################################################



###########################################################################
## Draw GO BP scatter plot: Up and Down DEG GOs together
goScatter = ggplot(data = goData) +
  geom_point(mapping = aes_string(x = "Term", 
                                  y = "richFactor", 
                                  size = "Count", 
                                  color = log10Field,
                                  fill = log10Field,
                                  shape = "up_or_down")) +
  scale_shape_manual(values =  c(24, 25), 
                     limits = c("up", "down"), 
                     labels = c("Up", "Down"),
                     name = "Differential expression\ntype") +
  scale_fill_gradientn(name = "-log10(p-value)",
                        colours = c("red", "green", "blue")) +
  scale_color_gradientn(colours = c("red", "green", "blue"), guide=FALSE) +
  scale_size_continuous(name = "Gene count", range = c(0, 10)) +
  scale_x_discrete(labels = wrap_format(80)) +
  guides(size = guide_legend(override.aes = list(shape = 17)),
         shape = guide_legend(override.aes = list(size = 6))) + 
  ggtitle(paste("Enriched GO terms (p-value <= 0.05) in", comparison, "comparison")) +
  ylab("Rich factor") + xlab("Biological Process GO Terms") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) + 
  coord_flip()


png(filename = paste0(outPrefix, "_GO_all_scatter.png", collapse = ""), width=12000, height=15000, res = 1000)
print(goScatter)
dev.off()

###########################################################################



###########################################################################
## Draw GO BP scatter plot: Up DEGs only
upGoData = dplyr::filter(goData, up_or_down == "up")

goTitle = paste("Enriched GO terms for up-regulated DEGs in", comparison, "comparison")
upGoScatter = topGO_scatterPlot(df = upGoData,
                                title = goTitle,
                                pvalCol = "log10_PValue",
                                termCol = "Term",
                                richCol = "richFactor",
                                geneCountCol = "Count")


png(filename = paste0(outPrefix, "_GO_up_scatter.png", collapse = ""), width = 4500, height = 5000, res = 400)
print(upGoScatter)
dev.off()

###########################################################################



###########################################################################
## Draw GO BP scatter plot: Down DEGs only
downGoData = dplyr::filter(goData, up_or_down == "down")

goTitle = paste("Enriched GO terms for down-regulated DEGs in", comparison, "comparison")

downGoScatter = topGO_scatterPlot(df = downGoData,
                                  title = goTitle,
                                  pvalCol = "log10_PValue",
                                  termCol = "Term",
                                  richCol = "richFactor",
                                  geneCountCol = "Count")


png(filename = paste0(outPrefix, "_GO_down_scatter.png", collapse = ""), width=4000, height=5000, res = 400)
print(downGoScatter)
dev.off()

###########################################################################



###########################################################################
## Draw KEGG pathway enrichment scatter plot: Up and Down DEG GOs together
pathData = data %>% dplyr::filter(Category == "KEGG_PATHWAY", UQ(as.name(pvalField)) <= 0.1) 
pathData$Term = gsub("dre\\d+:", "", pathData$Term, perl = T)

pathData = dplyr::arrange(pathData, up_or_down, !!as.name(pvalField))
pathData$up_or_down = factor(pathData$up_or_down)
pathData$Term = factor(pathData$Term, levels = unique(pathData$Term))


pathScatter = ggplot(data = pathData) +
  geom_point(mapping = aes_string(x = "Term", 
                                  y = "richFactor", 
                                  size = "Count", 
                                  color = log10Field,
                                  fill = log10Field,
                                  shape = "up_or_down")) +
  scale_shape_manual(values =  c(24, 25), 
                     limits = c("up", "down"), 
                     labels = c("Up", "Down"),
                     name = "Differential expression\ntype") +
  scale_fill_gradientn(name = "-log10(p-value)",
                       colours = c("red", "green", "blue")) +
  scale_color_gradientn(colours = c("red", "green", "blue"), guide=FALSE) +
  scale_size_continuous(name = "Gene count", range = c(0, 10)) +
  scale_x_discrete(labels = wrap_format(80)) +
  guides(size = guide_legend(override.aes = list(shape = 17)),
         shape = guide_legend(override.aes = list(size = 6))) + 
  ggtitle(paste("Enriched KEGG pathways in", comparison, "comparison")) +
  ylab("Rich factor") + xlab("KEGG pathways") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) + 
  coord_flip()


png(filename = paste0(outPrefix, "_KEGG_all_scatter.png", collapse = ""), width=8000, height=5000, res = 800)
print(pathScatter)
dev.off()

###########################################################################



###########################################################################
## Draw KEGG pathway enrichment scatter plot: Up DEG KEGG pathways
upKeggData = dplyr::filter(pathData, up_or_down == "up") %>%
  dplyr::arrange(up_or_down, !!as.name(pvalField))

goTitle = paste("Enriched KEGG pathways for up-regulated DEGs in", comparison, "comparison")
upKeggScatter = topGO_scatterPlot(df = upKeggData,
                                  title = goTitle,
                                  pvalCol = "log10_PValue",
                                  termCol = "Term",
                                  richCol = "richFactor",
                                  geneCountCol = "Count")


png(filename = paste0(outPrefix, "_KEGG_up_scatter.png", collapse = ""), width=4000, height=3000, res = 400)
print(upKeggScatter)
dev.off()

###########################################################################



###########################################################################
## Draw KEGG pathway enrichment scatter plot: Down DEG KEGG pathways
downKeggData = dplyr::filter(pathData, up_or_down == "down") %>%
  dplyr::arrange(up_or_down, !!as.name(pvalField))

goTitle = paste("Enriched KEGG pathways for down-regulated DEGs in", comparison, "comparison")
downKeggScatter = topGO_scatterPlot(df = downKeggData,
                                    title = goTitle,
                                    pvalCol = "log10_PValue",
                                    termCol = "Term",
                                    richCol = "richFactor",
                                    geneCountCol = "Count")

png(filename = paste0(outPrefix, "_KEGG_down_scatter.png", collapse = ""), width=4000, height=3000, res = 400)
print(downKeggScatter)
dev.off()


###########################################################################




###########################################################################
## Heatmap
# spread the table using dplyr package
goData = spread(data = goData, key = up_or_down, value = PValue) %>% dplyr::arrange(up, down)

rownames(goData) <- NULL
goData = column_to_rownames(goData, var = "Term") 

goDataMat = as.matrix(goData)

#Plot heatmap
ha1 = Heatmap(goDataMat, 
               col = colorRamp2(breaks = c(0, max(goDataMat, na.rm = T)), c("grey", "red"), space = "LAB"),
               column_title = "Enriched GO terms (p-value <= 0.05) in PG vs PV comparison",
               column_title_gp = gpar(fontface = "bold"),
               row_names_max_width = unit(15, "cm"), column_names_max_height = unit(6, "cm"), width = unit(8, "cm"),
               # row_dend_width = unit(3, "cm"), 
               cluster_columns = F, cluster_rows = F, 
               na_col = "grey", 
               heatmap_legend_param = list(title = "-LOG10(p-value)", color_bar = "continuous", title_position = 'topcenter', legend_direction = 'horizontal', legend_width = unit(5, "cm"), title_gp = gpar(fontsize = 12)))

draw(ha1, heatmap_legend_side = "bottom")


###########################################################################






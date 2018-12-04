
# library(ballgown)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(ggrepel)
library(data.table)
library(DESeq2)
library(tximport)

## This script uses ballgown to extract the FPKM matrix from the stringTie output
## It also perform PCA analysis using FPKM values
rm(list = ls())

path = "E:/Chris_UM/Analysis/28_RNAseq_mix_process/CL2018_RNAseq_01"

setwd(path)

file_sampleInfo = "sampleInfo_Ca.txt"
path_stringtie = "E:/Chris_UM/Analysis/28_RNAseq_mix_process/CL2018_RNAseq_01/stringTie"
## E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/stringTie/stringTie_CAlb
## E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/stringTie/stringTie_Caur


outPrefix = "CAlb_haploid_diff"

## samples to include. Should match with the values in 'condition' column in sampleInfo.txt
# compare = c("noDOX", "withDOX")

exptInfo = fread(file = file_sampleInfo, sep = "\t", header = T, stringsAsFactors = F)

# rownames(exptInfo) = exptInfo$sampleId

## set the factor levels.
## control levels should be first
# exptInfo$gt = factor(exptInfo$gt, levels = c("WT", "tet90"))
# exptInfo$treatment = factor(exptInfo$treatment, levels = c("YPD", "YPD_DOX", "YPD_GdA"))
exptInfo$condition = factor(exptInfo$condition, levels = c("SC5314", "GZ792", "GZ892"))


###########################################################################
## import the counts data using tximport and run DESeq2

design = ~ condition

filesStringtie = paste(path_stringtie, "/stringTie_", exptInfo$sampleId, "/t_data.ctab", sep = "")

names(filesStringtie) = exptInfo$sampleId

tmp = data.table::fread(file = filesStringtie[1], sep = "\t", header = T, stringsAsFactors = F)
tx2gene = tmp[, c("t_name", "gene_id")]

txi = tximport(files = filesStringtie, type = "stringtie", tx2gene = tx2gene)

ddsTxi = DESeqDataSetFromTximport(txi = txi, colData = exptInfo, design = design)

fpkmCounts = tibble::rownames_to_column(as.data.frame(fpkm(ddsTxi)), var = "geneId")

# assay(ddsTxi)
# colData(ddsTxi)
# rowData(ddsTxi)

## Run DESeq2
dds = DESeq(ddsTxi)

## raw counts
rawCounts = tibble::rownames_to_column(as.data.frame(counts(dds, normalized = FALSE)), var = "geneId")

fwrite(x = rawCounts, file = paste(outPrefix, "_rawCounts.tab", sep = ""),
       sep = "\t", row.names = F, col.names = T, quote = F)

## FPKM
fpkmCounts = tibble::rownames_to_column(as.data.frame(fpkm(dds)), var = "geneId")

fwrite(x = fpkmCounts, file = paste(outPrefix, "_FPKM.tab", sep = ""),
       sep = "\t", row.names = F, col.names = T, quote = F)


## r-log normalized counts
rld = rlog(dds)
rldCount = rownames_to_column(as.data.frame(assay(rld)), var = "geneID")

fwrite(x = rldCount, file = paste(outPrefix, "_rlogCounts.tab", sep = ""),
       sep = "\t", row.names = F, col.names = T, quote = F)


## PCA based on rld
plotPCA(rld, intgroup=c("condition"))

pcaData = plotPCA(rld, intgroup=c("condition"), returnData = TRUE)
percentVar = sprintf("%.2f", 100 * attr(pcaData, "percentVar"))

pltTitle = "Principal Component Analysis"
pointCol = base::structure(RColorBrewer::brewer.pal(n = length(unique(pcaData$condition)), name = "Set1"),
                           names = levels(pcaData$condition))


p1 = ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = condition), size=4) +
  geom_text_repel(mapping = aes(label = name), size = 3, point.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = pointCol) +
  ggtitle(pltTitle) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold", size = 15)
  )


png(filename = paste(outPrefix, "_PCA.png", sep = ""), width = 2000, height = 2000, res = 250)
p1
dev.off()


# ###########################################################################
# ## read the FPKM data from stringTie output using ballgown package
# 
# filesStringtie = paste(path_stringtie, "/stringTie_", exptInfo$sampleId, sep = "")
# 
# bg = ballgown(samples = filesStringtie, meas='all')
# 
# geneExpression = data.frame(gexpr(bg)) %>%
#   rownames_to_column(var = "geneId")
# 
# 
# ## rename the sample names
# renameCols = base::structure(gsub(pattern = "FPKM.stringTie_", replacement = "", x = names(geneExpression)),
#                              names = names(geneExpression))
# 
# 
# renameCols[colnames(geneExpression)]
# 
# colnames(geneExpression) = renameCols[colnames(geneExpression)]
# 
# 
# fwrite(x = geneExpression, file = paste(outPrefix, "FPKM_matrix.tab", sep = "_"),
#        sep = "\t", col.names = T, row.names = F, quote = F)
# 




###########################################################################
## PCA based on FPKM

## transform the data such that the genes are columns and each sample is a row
## also append the additional information for each sample using left_join()
df2 = data.table::dcast.data.table(data = data.table::melt.data.table(as.data.table(geneExpression), id.vars = "geneId"),
                                   formula = variable ~ geneId) %>% 
  dplyr::rename(sampleId = variable) %>% 
  dplyr::mutate(sampleId = as.character(sampleId)) %>% 
  dplyr::left_join(y = exptInfo, by = c("sampleId" = "sampleId")) %>% 
  dplyr::select(sampleId, gt, treatment, condition, replicate, dplyr::everything())


row.names(df2) = df2$sampleId

res.pca <- PCA(df2, graph = FALSE, scale.unit = TRUE, quali.sup = 1:5)

eig.val <- get_eigenvalue(res.pca)

## scree plot: variance by PC
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

## Graph of individuals
ind = get_pca_ind(res.pca)

fviz_pca_ind(res.pca,
             # col.ind = exprData$treatment,
             fill.ind = df2$condition,
             pointshape = 21,
             repel = TRUE,
             mean.point = FALSE,
             legend.title = "Study",
             pointsize = 3
)


## prepare the plot dataframe for ggplot
plotData = as.data.frame(ind$coord) %>%
  tibble::rownames_to_column(var = "sampleId") %>%
  dplyr::left_join(y = exptInfo, by = c("sampleId" = "sampleId"))

## set the factor levels
plotData$condition = factor(plotData$condition, levels = unique(exptInfo$condition))
plotData$treatment = factor(plotData$treatment, levels = levels(exptInfo$treatment))

pltTitle = "Principal Component Analysis"
pointCol = base::structure(RColorBrewer::brewer.pal(n = length(unique(plotData$condition)), name = "Set1"),
                           names = levels(plotData$condition))

pcaPlot = ggplot(data = plotData, mapping = aes(x = Dim.1, y = Dim.2, label = sampleId)) +
  geom_point(mapping = aes(color = condition),  size = 6) +
  geom_text_repel(size = 3, point.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = pointCol) +
  xlab( paste("PC1 (", sprintf("%.2f", eig.val[1, "variance.percent"]), "%)", sep = "") ) +
  ylab( paste("PC2 (", sprintf("%.2f", eig.val[2, "variance.percent"]), "%)", sep = "") ) +
  ggtitle(pltTitle) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"))



png(filename = paste(outPrefix, "_PCA.png", sep = ""), width = 4000, height = 3000, res = 420)
print(pcaPlot)
dev.off()


# ###########################################################################
# ## calculate mean FPKM
# 
# exptInfo = read.table(file = sampleInfoFile, header = T, sep = "\t", row.names = "sampleId")
# 
# ## IF NEEDED: modify the row names and column names
# # colnames(rawData) = sub("_WT", "_", colnames(rawData))
# # rownames(exptInfo) = sub("_WT", "_", rownames(exptInfo))
# 
# ## select only those sample rows which are part of current comparison
# designInfo = exptInfo
# # designInfo = droplevels(subset(exptInfo, condition %in% compare))
# 
# if(! all( rownames(designInfo) %in% colnames(geneExpression) ) ){
#   stop("Column names in FPKM matrix does not match with row names in experiment data")
# }
# 
# 
# ## get the sample names for each condition under comparison
# grp1 = sapply(rownames(designInfo[designInfo$condition %in% compare[1], ]), FUN = as.name, USE.NAMES = F, simplify = T)
# name1 = paste(compare[1], "_meanFPKM", sep = "")
# grp1Len = length(grp1)
# 
# grp2 = sapply(rownames(designInfo[designInfo$condition %in% compare[2], ]), FUN = as.name, USE.NAMES = F, simplify = T)
# name2 = paste(compare[2], "_meanFPKM", sep = "")
# grp2Len = length(grp2)
# 
# ## can add multiple groups
# 
# ## calculate the average FPKM
# fpkmData = geneExpression %>% 
#   dplyr::select(geneId, !!!c(grp1, grp2)) %>%
#   rowwise() %>%
#   mutate(!!name1 := sum(!!!grp1) / !!grp1Len,
#          !!name2 := sum(!!!grp2) / !!grp2Len)
# 
# 
# 
# ## write the data to file
# fwrite(x = fpkmData, file = "FPKM_matrix.tab", sep = "\t", quote = F, col.names = T)
# 
# 



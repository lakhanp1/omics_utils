library(DESeq2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(ggrepel)


## This script:
## 1. Runs the DESeq2 pipeline on count matrix (CSV) generated from StringTie
## 2. Writes the normalized count matrix and rlog transformed counts
## 3. Writes the differential expression results to Excel file
## 4. Creates volcano plot
##
## IMP: two values for log2foldChange are reported, before shrink and after shrinking.
## For selecting significant DEGs, the filter is applied on the log2foldChange values before shrinking
##

rm(list = ls())
source("E:/Chris_UM/Codes/RNA_Seq/DESeq2_functions.R")

path = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/CAuris_diff/"

## the denominator or WT in log2(fold_change) should be first
compare = c("WT_YPD", "WT_YPD_GdA")

outFilePrefix = paste("CAur_", paste(compare, collapse = "_vs_"), sep = "")
# outFilePrefix = "noDOX_withDOX_diff_retest"


path = paste(path, outFilePrefix, sep = "/")

if(!dir.exists(path)){
  dir.create(path = path)
}

setwd(path)


geneInfoFile = "E:/Chris_UM/Database/Candida_auris_B8441/Candida_auris_B8441.info.tab"
rawCountFile = "gene_count_matrix.csv"
sampleInfoFile = "sampleInfo.txt"




###########################################################################
##Prepare input data

p_cutoff <- 0.05
lfc_cut = 1
up_cut = lfc_cut
down_cut = lfc_cut * -1

rawData = fread(file = rawCountFile, header = T, sep = ",", stringsAsFactors = F, data.table = F)
rawData = column_to_rownames(rawData, var = "gene_id")

exptInfo = read.table(file = sampleInfoFile, header = T, sep = "\t", row.names = "sampleId", stringsAsFactors = F)


## IF NEEDED: modify the row names and column names
# colnames(rawData) = sub("_WT", "_", colnames(rawData))
# rownames(exptInfo) = sub("_WT", "_", rownames(exptInfo))


## select only those sample rows which are part of current comparison
designInfo = exptInfo
# designInfo = droplevels(subset(exptInfo, condition %in% compare))
countMat = rawData[, rownames(designInfo)]

## uncomment the 2 lines below if all conditions to be used in the analysis
# countMat <- rawData
# designInfo <- exptInfo
# subset(designInfo, condition %in% compare)

if(all(rownames(designInfo) %in% colnames(countMat))){
  countMat = countMat[, rownames(designInfo)]
} else{
  stop("Column names in count matrix does not match with row names in experiment data")
}

###########################################################################



###########################################################################
## Add gene symbol for each Ensembl ID

## Make sure to merge the geneIDs which are duplicated using group_by + summarize_all and remove the NA strings in data frame by replacing it with real <NA>
## geneSym <- read.table(file = ensToGeneId, sep = "\t", header = T, stringsAsFactors = F, na.strings = "") %>%
##   group_by(geneId) %>% summarize_all(.funs = funs(paste0(unique(.), collapse = ","))) %>% 
##   mutate_all(.funs = funs(ifelse(. == "NA", NA, .)))
##
## OR just take first row in case of duplicate rows 
geneSym = fread(file = geneInfoFile, sep = "\t", header = T, stringsAsFactors = F, na.strings = "") %>%
  distinct(geneId, .keep_all = T)

###########################################################################



###########################################################################
## run DESeq2 and extract the processed data
dds = DESeqDataSetFromMatrix(countData = countMat, colData = designInfo, design = ~ condition)

# colData(dds)

## Run DESeq2
dds = DESeq(dds)


## normalized counts matrix
normCounts = rownames_to_column(as.data.frame(counts(object = dds, normalized = T)), var = "geneID")

write.table(x = normCounts, file = paste0(c(outFilePrefix,"_normCounts.tab"), collapse = ""), row.names = F, col.names = T, sep = "\t", quote = F)


## r-log normalized counts
rld = rlog(dds)
rldCount = rownames_to_column(as.data.frame(assay(rld)), var = "geneID")

write.table(x = rldCount, file = paste0(c(outFilePrefix,"_rlogCounts.tab"), collapse = ""), row.names = F, col.names = T, sep = "\t", quote = F)

###########################################################################



###########################################################################
## PCA analysis
##NOTE: Typically, we recommend users to run samples from all groups together, and then use the contrast argument of the results function to extract comparisons of interest after fitting the model using DESeq. The model fit by DESeq estimates a single dispersion parameter for each gene, which defines how far we expect the observed count for a sample will be from the mean value from the model given its size factor and its condition group. However, for some datasets, exploratory data analysis (EDA) plots could reveal that one or more groups has much higher within-group variability than the others. This is case where, by comparing groups A and B separately - subsetting a DESeqDataSet to only samples from those two groups and then running DESeq on this subset - will be more sensitive than a model including all samples together.

plotPCA(rld, intgroup=c("condition"))

pcaData = plotPCA(rld, intgroup=c("condition"), returnData = TRUE)
percentVar = sprintf("%.2f", 100 * attr(pcaData, "percentVar"))

p1 = ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = condition), size=3) +
  geom_text_repel(mapping = aes(label = name),
                  min.segment.length = unit(1, "lines")) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle(paste("PCA plot of all samples in comparison", paste0(compare, collapse = " vs "), sep = " ")) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold", size = 15)
  )


png(filename = "PCA_plot.png", width = 2000, height = 2000, res = 250)
p1
dev.off()


###########################################################################



###########################################################################
## result for the comparison of interest
## "Independent Filtering": For weakly expressed genes, we have no chance of seeing differential expression, because the low read counts suffer from such high Poisson noise that any biological effect is drowned in the uncertainties from the sampling at a low rate. Genes with very low mean count have little or no power, and are best excluded from testing. At first sight, there may seem to be little benefit in filtering out these genes. After all, the test found them to be non-significant anyway. However, these genes have an influence on the multiple testing adjustment, whose performance improves if such genes are removed. By removing the low count genes from the input to the FDR procedure, we can find more genes to be significant among those that we keep, and so improved the power of our test. This approach is known as independent filtering. REF: https://www.bioconductor.org/help/workflows/rnaseqGene/#independent-filtering

res = results(dds, contrast = c("condition", rev(compare)), alpha = 0.05)
summary(res)

## get shrunken log2 fold changes
resShrink = lfcShrink(dds, contrast = c("condition", rev(compare)), res=res)
summary(resShrink)
# mcols(resShrink, use.names=TRUE)


###########################################################################
## MA plot
png(filename = paste0(outFilePrefix, "_MA.png", collapse = ""), width = 3000, height = 3000, res = 300)
plotMA(res, ylim=c(-4,4), main = "MA plot")
dev.off()


png(filename = paste0(outFilePrefix, "_MA_shrunk.png", collapse = ""), width = 3000, height = 3000, res = 300)
plotMA(resShrink, ylim=c(-4,4), main = "MA plot with shrunken log2 fold changes")
dev.off()

# plotDispEsts(dds)

###########################################################################

## store un-shrunk LFC data
resDf = rownames_to_column(as.data.frame(res), var = "geneID")

## store shrunk LFC data
resShrinkDf = rownames_to_column(as.data.frame(resShrink), var = "geneID") %>%
  dplyr::select(geneID, shrinkLog2FC = log2FoldChange, sLfcSE = lfcSE)

resultTable = dplyr::left_join(x = resDf, y = resShrinkDf, by = c("geneID")) %>%
  dplyr::select(geneID, baseMean, log2FoldChange, lfcSE, shrinkLog2FC, sLfcSE, everything())


## Extract significant genes and write to Excel file
## get the sample names for each condition under comparison
grp1 = sapply(rownames(designInfo[designInfo$condition %in% compare[1], ]), FUN = as.name, USE.NAMES = F, simplify = T)
name1 = paste(compare[1], "_meanCount", sep = "")
grp1Len = length(grp1)

grp2 = sapply(rownames(designInfo[designInfo$condition %in% compare[2], ]), FUN = as.name, USE.NAMES = F, simplify = T)
name2 = paste(compare[2], "_meanCount", sep = "")
grp2Len = length(grp2)


resNormCounts = normCounts %>% 
  dplyr::select(geneID, !!!c(grp1, grp2)) %>%
  rowwise() %>%
  mutate(!!name1 := sum(!!!grp1) / !!grp1Len,
         !!name2 := sum(!!!grp2) / !!grp2Len)
  

## validate
# as.data.frame(rowMeans(head(resNormCounts)[,grp1]), col.names = c("grp1_mean"))


selCols = colnames(geneSym)[colnames(geneSym) != "geneId"]

diffData = resultTable %>% 
  left_join(y = resNormCounts, by = "geneID") %>% 
  left_join(y = geneSym, by = c("geneID" = "geneId")) %>%
  dplyr::select(geneID, !!!selCols, ends_with("_meanCount"), log2FoldChange, shrinkLog2FC, pvalue, padj) %>%
  mutate(
    diff_l2fc = case_when(
      padj < p_cutoff & log2FoldChange >= up_cut ~ "up",
      padj < p_cutoff & log2FoldChange <= down_cut ~ "down",
      TRUE ~ "noDEG"
    ),
    diff_shrink_l2fc = case_when(
      padj < p_cutoff & shrinkLog2FC >= up_cut ~ "up",
      padj < p_cutoff & shrinkLog2FC <= down_cut ~ "down",
      TRUE ~ "noDEG"
    )
  )


head(diffData)


significant = filter(diffData, padj < p_cutoff, log2FoldChange >= up_cut | log2FoldChange <= down_cut)

significant_up = filter(diffData, padj < p_cutoff, log2FoldChange >= up_cut)
significant_down = filter(diffData, padj < p_cutoff, log2FoldChange <= down_cut)


write.table(x = diffData, file = paste0(outFilePrefix, "_DEG_all.txt", collapse = ""), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(x = significant_up, file = paste0(outFilePrefix, "_DEG_up.txt", collapse = ""), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(x = significant_down, file = paste0(outFilePrefix, "_DEG_down.txt", collapse = ""), col.names = T, row.names = F, sep = "\t", quote = F)

###########################################################################








###########################################################################
## plot volcano plot
# diffData = fread(file = paste0(outFilePrefix, "_DEG_all.txt", collapse = ""), sep = "\t", header = T, stringsAsFactors = F, data.table = F)

markGenes = c()

# tmpDf = filter(diffData, geneName %in% markGenes)

plotTitle = paste("Volcano plot:", compare[1], "vs", compare[2], sep = " ")
p2 = volcanoPlot(df = diffData, 
                 namePrefix = outFilePrefix, 
                 title = plotTitle, 
                 fdr_col = "padj", 
                 lfc_col = "log2FoldChange", 
                 geneOfInterest = markGenes,
                 ylimit = 100)

png(filename = paste0(outFilePrefix, "_volcano.png", collapse = ""), width = 3000, height = 3000, res = 230)
p2
dev.off()


###########################################################################









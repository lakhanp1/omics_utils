library(DESeq2)
library(tidyverse)
library(data.table)
library(ggrepel)
library(tximport)
library(here)
library(org.Hsapiens.eg.db)


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

source("E:/Chris_UM/GitHub/omics_util/RNAseq_scripts/02_DESeq2_functions.R")

analysisName <- "MHCC97L_AB_vs_B"

## the denominator or WT in log2(fold_change) should be second
compare <- c("AB", "B")

file_sampleInfo <- here::here("data", "sample_info.txt")


outDir <- here::here("analysis", "02_DESeq2_diff", analysisName)
outPrefix <- paste(outDir, analysisName, sep = "/")
orgDb <- org.Hsapiens.eg.db


if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


FDR_cut <- 0.05
lfc_cut <- 0.585
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

###########################################################################
## set levels for experiment design information

exptInfo <- read.table(file = file_sampleInfo, header = T, sep = "\t", stringsAsFactors = F)

## set the reference levels
exptInfo$genotype <- factor(exptInfo$genotype, levels = c("MHCC97L"))
# "WT", "A", "B", "AB"
exptInfo$condition <- factor(
  x = exptInfo$condition,
  levels = c("B", "A", "WT", "AB")
)

rownames(exptInfo) <- exptInfo$sampleId


## select only those sample rows which are part of current comparison
# exptInfo <- droplevels(subset(exptInfo, condition %in% compare))

design <- ~ condition

###########################################################################
## import counts data: either by tximport or as raw count matrix

## import the counts data using tximport and run DESeq2
path_stringtie <- here::here("data", "stringTie")
filesStringtie <- paste(path_stringtie, "/stringTie_", exptInfo$sampleId, "/t_data.ctab", sep = "")
names(filesStringtie) <- exptInfo$sampleId

tmp <- data.table::fread(file = filesStringtie[1], sep = "\t", header = T, stringsAsFactors = F)
tx2gene <- tmp[, c("t_name", "gene_id")]

txi <- tximport(files = filesStringtie, type = "stringtie", tx2gene = tx2gene,
                readLength = 100)

ddsTxi <- DESeqDataSetFromTximport(txi = txi, colData = exptInfo, design = design)
# assay(ddsTxi)
# colData(ddsTxi)
# rowData(ddsTxi)

## Run DESeq2
dds <- DESeq(ddsTxi)


# ## import raw counts data and run DESeq2
# file_rawCounts <- here::here("RNAseq_data", "MatrixCountsPerGeneBySample.txt")
# 
# countsDf <- suppressMessages(readr::read_tsv(file = file_rawCounts, col_names = T)) %>%
#   as.data.frame()
# rownames(countsDf) <- countsDf$geneId
# countsDf$geneId <- NULL
# 
# 
# if(all(rownames(exptInfo) %in% colnames(countsDf))){
#   countsDf <- countsDf[, rownames(exptInfo)]
# } else{
#   stop("Column names in count matrix does not match with row names in experiment data")
# }
# 
# ## select only those sample rows which are part of current comparison
# # exptInfo <- droplevels(subset(exptInfo, condition %in% compare))
# # countsDf <- countsDf[, rownames(exptInfo)]
# 
# ## run DESeq2 and extract the processed data
# ddsCount <- DESeqDataSetFromMatrix(countData = countsDf, colData = exptInfo, design = design)
# 
# ## Run DESeq2
# dds <- DESeq(ddsCount)

###########################################################################
## Add gene symbol for each Ensembl ID

## Make sure to merge the geneIds which are duplicated using group_by + summarize_all and remove the NA strings in data frame by replacing it with real <NA>
## geneInfo <- read.table(file = ensToGeneId, sep = "\t", header = T, stringsAsFactors = F, na.strings = "") %>%
##   group_by(geneId) %>% summarize_all(.funs = funs(paste0(unique(.), collapse = ","))) %>% 
##   mutate_all(.funs = funs(ifelse(. == "NA", NA, .)))
##
## OR just take first row in case of duplicate rows 
# geneInfo <- fread(file = geneInfoFile, sep = "\t", header = T, stringsAsFactors = F, na.strings = "") %>%
#   distinct(geneId, .keep_all = T)


## use org.db
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(x = orgDb, keytype = "ENSEMBL_VERSION"),
                                  columns = c("GENE_NAME", "DESCRIPTION"),
                                  keytype = "ENSEMBL_VERSION") %>% 
  dplyr::rename(geneId = ENSEMBL_VERSION)


###########################################################################

## raw counts
rawCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = FALSE)), var = "geneId")

readr::write_tsv(x = rawCounts, path = paste(outPrefix, ".rawCounts.tab", sep = ""))

## FPKM
# fpkmCounts <- tibble::rownames_to_column(as.data.frame(fpkm(dds)), var = "geneId")
# 
# fwrite(x = fpkmCounts, file = paste(outPrefix, "_FPKM.tab", sep = ""),
#        sep = "\t", row.names = F, col.names = T, quote = F)


## normalized counts matrix
normCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = TRUE)), var = "geneId")

readr::write_tsv(x = normCounts, path = paste0(c(outPrefix,".normCounts.tab"), collapse = ""))


## r-log normalized counts
rld <- rlog(dds, blind = FALSE)
rldCount <- rownames_to_column(as.data.frame(assay(rld)), var = "geneId")

readr::write_tsv(x = rldCount, path = paste(outPrefix, ".rlogCounts.tab", sep = ""))


###########################################################################
## PCA analysis
##NOTE: Typically, we recommend users to run samples from all groups together, and then use the contrast argument of the results function to extract comparisons of interest after fitting the model using DESeq. The model fit by DESeq estimates a single dispersion parameter for each gene, which defines how far we expect the observed count for a sample will be from the mean value from the model given its size factor and its condition group. However, for some datasets, exploratory data analysis (EDA) plots could reveal that one or more groups has much higher within-group variability than the others. This is case where, by comparing groups A and B separately - subsetting a DESeqDataSet to only samples from those two groups and then running DESeq on this subset - will be more sensitive than a model including all samples together.


plotPCA(rld, intgroup=c("condition"))

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData = TRUE)
percentVar <- sprintf("%.2f", 100 * attr(pcaData, "percentVar"))

pltTitle <- paste("Principal Component Analysis:", compare[1], "vs", compare[2])
pointCol <- base::structure(RColorBrewer::brewer.pal(n = length(unique(pcaData$condition)), name = "Set1"),
                            names = levels(pcaData$condition))


pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
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
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold", size = 15)
  )


png(filename = paste(outPrefix, ".PCA.png", sep = ""), width = 3000, height = 3000, res = 300)
pcaPlot
dev.off()

###########################################################################



###########################################################################
## result for the comparison of interest
## "Independent Filtering": For weakly expressed genes, we have no chance of seeing differential expression, because the low read counts suffer from such high Poisson noise that any biological effect is drowned in the uncertainties from the sampling at a low rate. Genes with very low mean count have little or no power, and are best excluded from testing. At first sight, there may seem to be little benefit in filtering out these genes. After all, the test found them to be non-significant anyway. However, these genes have an influence on the multiple testing adjustment, whose performance improves if such genes are removed. By removing the low count genes from the input to the FDR procedure, we can find more genes to be significant among those that we keep, and so improved the power of our test. This approach is known as independent filtering. REF: https://www.bioconductor.org/help/workflows/rnaseqGene/#independent-filtering

resultsNames(dds)
# coefName <- "condition_5A9_Control_vs_CEA17_Control"
coefName <- paste("condition_", compare[1], "_vs_", compare[2], sep = "")

if(!any(coefName %in% resultsNames(dds))){
  stop("Wrong coefficient name")
}

res <- results(dds, name = coefName, alpha = 0.05)
summary(res)

## get shrunken log2 fold changes
resShrink <- lfcShrink(dds, coef = coefName,
                       res=res, type="apeglm")

summary(resShrink)
mcols(resShrink, use.names=TRUE)


fcDensity <- hist(x = pmin(pmax(res$log2FoldChange, -5), 5),
                  breaks = 50,
                  main = "log2(fold-change) frequency distribution")

hist(x = res$pvalue, breaks = 100, main = "p-value distribution")
hist(x = res$pvalue, breaks = c(seq(0, 0.1, length.out = 20), seq(0.11,1, by = 0.01)),
     main = "p-value distribution")

hist(x = res$padj, breaks = 100, main = "q-value distribution")

lfcFreq <- ggplot(data = as.data.frame(res)) +
  geom_histogram(
    mapping = aes(x = pmin(pmax(log2FoldChange, -5), 5)),
    fill = "red", bins = 50) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  scale_x_continuous(limits = c(-5, 5), expand = expand_scale(mult = 0.01)) +
  scale_y_continuous(expand = expand_scale(mult = 0.01)) +
  labs(title = paste("log2(fold change) frequency distribution:", compare[1], "vs", compare[2]),
       x = "log2(fold change)", y = "Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(face = "bold", size = 15),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = "none"
  )


pqDensity <- ggplot(data = as.data.frame(res)) +
  geom_histogram(mapping = aes(x = pvalue, y= ..density.. , fill = "pvalue"),
                 bins = 100, alpha = 0.5) +
  geom_histogram(mapping = aes(x = padj, y= ..density.. , fill = "padj"),
                 bins = 100, alpha = 0.5) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  scale_fill_manual(
    name = NULL,
    values = c("pvalue" = "#999999", "padj" = "#E69F00"), 
    breaks = c("pvalue", "padj"),
    labels = c("p-value", "q-value")
  ) +
  scale_x_continuous(expand = expand_scale(mult = 0.01)) +
  scale_y_continuous(expand = expand_scale(mult = 0.01)) +
  labs(title = paste("p-value and q-value density distribution:", compare[1], "vs", compare[2]),
       x = "p-value or q-value", y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.2,"cm"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.title = element_text(face = "bold", size = 15),
        legend.position = c(0.5,0.8)
  )


###########################################################################
## MA plot
# png(filename = paste(outPrefix, ".MA.png", sep = ""), width = 2000, height = 3000, res = 250)
# pdf(file = paste(outPrefix, ".MA.pdf", sep = ""), width = 8, height = 10, onefile = F)
# op <- par(mfrow = c(2, 1))

plotMA(res, ylim=c(-4,4), main = "MA plot with unshrunken LFC")
plotMA(resShrink, ylim=c(-4,4), main = "MA plot with shrunken log2 fold changes")


# ## MA plot with unshrunken LFC
# geneplotter::plotMA(object = data.frame(
#   m = res$baseMean,
#   a = res$log2FoldChange,
#   c = ifelse(abs(res$log2FoldChange) >= lfc_cut & res$padj <= FDR_cut, TRUE, FALSE)),
#   ylim=c(-4,4), 
#   main = "MA plot with unshrunken LFC"
# )
# 
# ## MA plot with apeglm shrunken LFC
# geneplotter::plotMA(object = data.frame(
#   m = resShrink$baseMean,
#   a = resShrink$log2FoldChange,
#   c = ifelse(abs(resShrink$log2FoldChange) >= lfc_cut & resShrink$padj <= FDR_cut, TRUE, FALSE)),
#   ylim=c(-4,4), 
#   main = "MA plot with apeglm shrunken LFC"
# )

# par(op)
# dev.off()


###########################################################################

## store un-shrunk LFC data
resDf <- rownames_to_column(as.data.frame(res), var = "geneId")

## store shrunk LFC data
resShrinkDf <- rownames_to_column(as.data.frame(resShrink), var = "geneId") %>% 
  dplyr::mutate(contrast = coefName)

resultTable <- dplyr::left_join(
  x = resDf,
  y = dplyr::select(resShrinkDf, geneId, contrast, shrinkLog2FC = log2FoldChange, sLfcSE = lfcSE),
  by = c("geneId")) %>%
  dplyr::select(geneId, contrast, baseMean, log2FoldChange, lfcSE, shrinkLog2FC, sLfcSE, everything())


## Extract significant genes and write to Excel file
## get the sample names for each condition under comparison
grp1 <- sapply(rownames(exptInfo[exptInfo$condition %in% compare[1], ]), FUN = as.name, USE.NAMES = F, simplify = T)
name1 <- paste(compare[1], "_meanCount", sep = "")
grp1Len <- length(grp1)

grp2 <- sapply(rownames(exptInfo[exptInfo$condition %in% compare[2], ]), FUN = as.name, USE.NAMES = F, simplify = T)
name2 <- paste(compare[2], "_meanCount", sep = "")
grp2Len <- length(grp2)


resNormCounts <- normCounts %>% 
  dplyr::select(geneId, !!!c(grp1, grp2)) %>%
  rowwise() %>%
  mutate(!!name1 := sum(!!!grp1) / !!grp1Len,
         !!name2 := sum(!!!grp2) / !!grp2Len)


## validate
# as.data.frame(rowMeans(head(resNormCounts)[,grp1]), col.names = c("grp1_mean"))


selCols <- colnames(geneInfo)[colnames(geneInfo) != "geneId"]

diffData <- resultTable %>% 
  left_join(y = resNormCounts, by = "geneId") %>% 
  left_join(y = geneInfo, by = c("geneId" = "geneId")) %>%
  mutate(
    diff_l2fc = dplyr::case_when(
      padj < FDR_cut & log2FoldChange >= up_cut ~ "up",
      padj < FDR_cut & log2FoldChange <= down_cut ~ "down",
      TRUE ~ "noDEG"
    ),
    diff_shrink_l2fc = dplyr::case_when(
      padj < FDR_cut & shrinkLog2FC >= up_cut ~ "up",
      padj < FDR_cut & shrinkLog2FC <= down_cut ~ "down",
      TRUE ~ "noDEG"
    )
  ) %>% 
  dplyr::select(geneId, contrast, ends_with("_meanCount"), log2FoldChange, shrinkLog2FC,
                pvalue, padj, diff_l2fc, diff_shrink_l2fc, !!!selCols)


head(diffData)


significant <- filter(diffData, padj < FDR_cut, log2FoldChange >= up_cut | log2FoldChange <= down_cut)

significant_up <- filter(diffData, padj < FDR_cut, log2FoldChange >= up_cut)
significant_down <- filter(diffData, padj < FDR_cut, log2FoldChange <= down_cut)


readr::write_tsv(x = resShrinkDf, path = paste(outPrefix, ".DESeq2.tab", sep = ""))
readr::write_tsv(x = diffData, path = paste(outPrefix, ".DEG_all.txt", sep = ""))


###########################################################################
## plot volcano plot
# diffData <- fread(file = paste0(outPrefix, "_DEG_all.txt", collapse = ""), sep = "\t", header = T, stringsAsFactors = F, data.table = F)

markGenes <- c()

# tmpDf = filter(diffData, geneName %in% markGenes)

plotTitle <- paste("Volcano plot:", compare[1], "vs", compare[2], sep = " ")
p2 <- volcano_plot(df = diffData, 
                   title = plotTitle, 
                   fdr_col = "padj", 
                   lfc_col = "log2FoldChange",
                   fdr_cut = FDR_cut, lfc_cut = lfc_cut,
                   geneOfInterest = markGenes,
                   ylimit = 100, xlimit = c(-7, 7))

png(filename = paste(outPrefix, ".volcano.png", sep = ""), width = 2500, height = 3000, res = 280)
plot(p2$plot)
dev.off()

###########################################################################

# plot all data in single PDF file

pdf(file = paste(outPrefix, ".summary_plots.pdf", sep = ""), width = 10, height = 10, onefile = TRUE)

## PCA
plot(pcaPlot)

## p-value distribution plots
plot(pqDensity)
plot(lfcFreq)

## MA plots
op <- par(mfrow = c(2, 1))
plotMA(
  object = res,
  ylim=c(-4,4),
  main = paste("MA plot with unshrunken LFC:", compare[1], "vs", compare[2])
)

plotMA(
  object = resShrink,
  ylim=c(-4,4),
  main = paste("MA plot with shrunken LFC:", compare[1], "vs", compare[2])
)

par(op)

## volcano plot
plot(p2$plot)

dev.off()

###########################################################################



library(DESeq2)
library(tidyverse)
library(data.table)
library(ggrepel)
library(tximport)
library(org.AFumigatus293.eg.db)
library(here)

## This script:
## read the stringtie output using tximport or as raw counts
## perform differential gene expression analysis with interaction terms using DESeq2

rm(list = ls())
source("E:/Chris_UM/GitHub/omics_util/RNAseq_scripts/DESeq2_functions.R")

analysisName <- "creE_del_AA_effect"

outDir <- here::here("RNAseq_data", analysisName)
file_sampleInfo <- here::here("RNAseq_data", "sampleInfo.txt")

outPrefix <- paste(outDir, analysisName, sep = "/")

file_geneInfo <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r09_geneInfo.tab"
# "E:/Chris_UM/Database/Candida_auris_B8441/Candida_auris_B8441.info.tab"
# "E:/Chris_UM/Database/C_albicans/SC5314_A21/C_albicans_SC5314_A21.info.tab"
# "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r09_geneInfo.tab"

orgDb <- org.AFumigatus293.eg.db
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
exptInfo$genotype <- factor(exptInfo$genotype, levels = c("CEA17", "5A9"))
exptInfo$treatment <- factor(exptInfo$treatment, levels = c("C", "AA"))
exptInfo$condition <- factor(exptInfo$condition, levels = c("CEA17_C", "5A9_C", "CEA17_AA", "5A9_AA"))

rownames(exptInfo) <- exptInfo$sampleId

design <- ~ genotype + treatment + genotype:treatment

###########################################################################
## Add gene information

# geneSym = data.table::fread(file = file_geneInfo, sep = "\t", header = T, stringsAsFactors = F) %>%
#   distinct(geneId, .keep_all = T)

###########################################################################
## import counts data: either by tximport or as raw count matrix

# ## import the counts data using tximport and run DESeq2
# path_stringtie <- here::here("data", "stringTie", "stringTie_CAur")
# filesStringtie <- paste(path_stringtie, "/stringTie_", exptInfo$sampleId, "/t_data.ctab", sep = "")
# names(filesStringtie) <- exptInfo$sampleId
# 
# tmp <- data.table::fread(file = filesStringtie[1], sep = "\t", header = T, stringsAsFactors = F)
# tx2gene <- tmp[, c("t_name", "gene_id")]
# 
# txi <- tximport(files = filesStringtie, type = "stringtie", tx2gene = tx2gene)
# 
# ddsTxi <- DESeqDataSetFromTximport(txi = txi, colData = exptInfo, design = design)
# assay(ddsTxi)
# colData(ddsTxi)
# rowData(ddsTxi)
# 
# ## Run DESeq2
# dds <- DESeq(ddsTxi)


## import raw counts data and run DESeq2
file_rawCounts <- here::here("RNAseq_data", "MatrixCountsPerGeneBySample.txt")

countsDf <- suppressMessages(readr::read_tsv(file = file_rawCounts, col_names = T)) %>%
  as.data.frame()
rownames(countsDf) <- countsDf$geneId
countsDf$geneId <- NULL

## select only those sample rows which are part of current comparison
# exptInfo <- droplevels(subset(exptInfo, condition %in% compare))
countsDf <- countsDf[, rownames(exptInfo)]

## run DESeq2 and extract the processed data
ddsCount <- DESeqDataSetFromMatrix(countData = countsDf, colData = exptInfo, design = design)

## Run DESeq2
dds <- DESeq(ddsCount)


###########################################################################
## raw counts
rawCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = FALSE)), var = "geneId")

readr::write_tsv(x = rawCounts, path = paste(outPrefix, ".rawCounts.tab", sep = ""))

## FPKM
# fpkmCounts <- tibble::rownames_to_column(as.data.frame(fpkm(dds)), var = "geneId")
# 
# fwrite(x = fpkmCounts, file = paste(outPrefix, ".FPKM.tab", sep = ""),
#        sep = "\t", row.names = F, col.names = T, quote = F)


## normalized counts matrix
normCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = TRUE)), var = "geneId")

readr::write_tsv(x = normCounts, path = paste0(c(outPrefix,".normCounts.tab"), collapse = ""))


## r-log normalized counts
rld <- rlog(dds)
rldCount <- rownames_to_column(as.data.frame(assay(rld)), var = "geneId")

readr::write_tsv(x = rldCount, path = paste(outPrefix, ".rlogCounts.tab", sep = ""))


###########################################################################
## PCA analysis
##NOTE: Typically, we recommend users to run samples from all groups together, and then use the contrast argument of the results function to extract comparisons of interest after fitting the model using DESeq. The model fit by DESeq estimates a single dispersion parameter for each gene, which defines how far we expect the observed count for a sample will be from the mean value from the model given its size factor and its condition group. However, for some datasets, exploratory data analysis (EDA) plots could reveal that one or more groups has much higher within-group variability than the others. This is case where, by comparing groups A and B separately - subsetting a DESeqDataSet to only samples from those two groups and then running DESeq on this subset - will be more sensitive than a model including all samples together.

plotPCA(rld, intgroup=c("condition"))

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData = TRUE)
percentVar <- sprintf("%.2f", 100 * attr(pcaData, "percentVar"))

pltTitle <- "Principal Component Analysis"
pointCol <- base::structure(RColorBrewer::brewer.pal(n = length(unique(pcaData$condition)), name = "Set1"),
                            names = levels(pcaData$condition))


p1 <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
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


pdf(file = paste(outPrefix, ".PCA.pdf", sep = ""), width = 10, height = 10)
p1
dev.off()


###########################################################################
## extract results for interaction terms

resultsNames(dds)

## the condition effect B_vs_A for genotype I (the main effect)
res_BA_gt1 <- results(dds, name = "treatment_AA_vs_C")

## another way to extract condition effect for genotype I using contrast argument
# res_BA_gt1_2 <- results(dds, contrast = c("treatment", "AA", "Control"), lfcThreshold = lfc_cut)

summary(res_BA_gt1)
mcols(res_BA_gt1)
plotMA(res_BA_gt1, ylim=c(-4,4), main = "MA plot with unshrunken LFC")

resLfcShrink_BA_gt1 <- lfcShrink(dds, coef = "treatment_AA_vs_C",
                                 res = res_BA_gt1,
                                 type="apeglm")

plotMA(resLfcShrink_BA_gt1, ylim=c(-4,4), main = "MA plot with shrunken LFC")

readr::write_tsv(x = rownames_to_column(as.data.frame(resLfcShrink_BA_gt1), var = "geneId"),
                 path = paste(outPrefix, ".AA_vs_C.CEA17.DESeq2.tab", sep = ""))


## the condition effect B_vs_A for genotype II
## remember that the contrast has to be a list with first element as a vector of genotype I effect and interaction term
res_BA_gt2 <- results(dds,
                      contrast = list(c("treatment_AA_vs_C", "genotype5A9.treatmentAA")))

summary(res_BA_gt2)
mcols(res_BA_gt2)
plotMA(res_BA_gt2, ylim=c(-4,4), main = "MA plot with unshrunken LFC")

## contrast cannot be used with apeglm. so use ashr for shrinking LFC
resLfcShrin_BA_gt2 <- lfcShrink(dds, contrast = list(c("treatment_AA_vs_C", "genotype5A9.treatmentAA")),
                                res = res_BA_gt2, type = "ashr")

plotMA(resLfcShrin_BA_gt2, ylim=c(-4,4), main = "MA plot with shrunken LFC")


readr::write_tsv(x = rownames_to_column(as.data.frame(resLfcShrin_BA_gt2), var = "geneId"),
                 path = paste(outPrefix, ".AA_vs_C.5A9.DESeq2.tab", sep = ""))

## the interaction term for condition effect B_vs_A in genotype II vs genotype I.
## the interaction term, answering: is the condition effect *different* across genotypes
inter_gt21_BA <- "genotype5A9.treatmentAA"

resInter_gt21_BA <- results(dds, name = inter_gt21_BA)

summary(resInter_gt21_BA)
mcols(resInter_gt21_BA)


resShrink_inter_gt21_BA <- lfcShrink(dds, coef = inter_gt21_BA,
                                     res = resInter_gt21_BA,
                                     type="apeglm")

summary(resShrink_inter_gt21_BA)
mcols(resShrink_inter_gt21_BA)




###########################################################################
## MA plot
# png(filename = paste(outPrefix, ".MA.png", sep = ""), width = 2000, height = 3000, res = 250)
pdf(file = paste(outPrefix, ".MA.pdf", sep = ""), width = 8, height = 8, onefile = F)
op <- par(mfrow = c(2, 1))

plotMA(resInter_gt21_BA, ylim=c(-4,4), main = "MA plot with unshrunken LFC")
plotMA(resShrink_inter_gt21_BA, ylim=c(-4,4), main = "MA plot with apeglm shrunken LFC")

# ## MA plot with unshrunken LFC
# geneplotter::plotMA(object = data.frame(
#   m = resInter_gt21_BA$baseMean,
#   a = resInter_gt21_BA$log2FoldChange,
#   c = ifelse(abs(resInter_gt21_BA$log2FoldChange) >= lfc_cut & resInter_gt21_BA$padj <= FDR_cut, TRUE, FALSE)),
#   ylim=c(-4,4), 
#   main = "MA plot with unshrunken LFC"
# )
# 
# ## MA plot with apeglm shrunken LFC
# geneplotter::plotMA(object = data.frame(
#   m = resShrink_inter_gt21_BA$baseMean,
#   a = resShrink_inter_gt21_BA$log2FoldChange,
#   c = ifelse(abs(resShrink_inter_gt21_BA$log2FoldChange) >= lfc_cut & resShrink_inter_gt21_BA$padj <= FDR_cut, TRUE, FALSE)),
#   ylim=c(-4,4), 
#   main = "MA plot with apeglm shrunken LFC"
# )

par(op)
dev.off()




###########################################################################
## generate the results DF

## store un-shrunk LFC data
resDf <- rownames_to_column(as.data.frame(resInter_gt21_BA), var = "geneId")
head(resDf)

## store shrunk LFC data
resShrinkDf <- rownames_to_column(as.data.frame(resShrink_inter_gt21_BA), var = "geneId") 
head(resShrinkDf)

resultTable <- dplyr::left_join(
  x = resDf,
  y = dplyr::select(resShrinkDf, geneId, shrinkLog2FC = log2FoldChange, sLfcSE = lfcSE),
  by = c("geneId")) %>%
  dplyr::select(geneId, baseMean, log2FoldChange, lfcSE, shrinkLog2FC, sLfcSE, everything())

# resultTable <- resShrinkDf
head(resultTable)



## function to get mean count formula as list
formula_list <- function(x){
  ids = samples = sapply(x$sampleId, FUN = as.name, USE.NAMES = F, simplify = T)
  len = nrow(x)
  name = as.name(paste(x$condition[1], "_meanCount", sep = ""))
  
  form = quos(!! name := sum(!!! ids) / !! len)
  
  return(form)
}


## generate formula list
af <- exptInfo %>% 
  # dplyr::mutate_if(.predicate = is.factor, .funs = as.character) %>% 
  dplyr::group_by(condition) %>% 
  dplyr::do(form = formula_list(.)) %>% 
  dplyr::mutate(name = paste(condition, "_meanCount", sep = ""))

avgFormula <- unlist(af$form)

## apply formula on each row to count mean count
countsDf <- dplyr::rowwise(normCounts) %>%
  dplyr::mutate(!!! avgFormula)

# write.table(head(countsDf), file = "clipboard", sep = "\t", row.names = F, quote = F)

geneSym <- AnnotationDbi::select(x = orgDb,
                                 keys = keys(orgDb, keytype = "GID"),
                                 columns = c("DESCRIPTION"), keytype = "GID") %>% 
  dplyr::rename(geneId = GID)

selCols <- colnames(geneSym)[colnames(geneSym) != "geneId"]

diffData <- resultTable %>% 
  left_join(y = countsDf, by = "geneId") %>% 
  left_join(y = geneSym, by = c("geneId" = "geneId")) %>%
  mutate(
    diff_l2fc = case_when(
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
  dplyr::select(geneId, ends_with("_meanCount"), log2FoldChange, shrinkLog2FC,
                pvalue, padj, diff_l2fc, diff_shrink_l2fc, !!!selCols)

head(diffData)

significant <- filter(diffData, padj < FDR_cut, log2FoldChange >= up_cut | log2FoldChange <= down_cut)

significant_up <- filter(diffData, padj < FDR_cut, log2FoldChange >= up_cut)
significant_down <- filter(diffData, padj < FDR_cut, log2FoldChange <= down_cut)

readr::write_tsv(x = resultTable, path = paste(outPrefix, ".DESeq2.tab", sep = ""))
readr::write_tsv(x = diffData, path = paste(outPrefix, ".DEG_all.txt", sep = ""))

# readr::write_tsv(x = significant_up, path = paste(outPrefix, "_DEG_up.txt", sep = ""))
# readr::write_tsv(x = significant_down, path = paste(outPrefix, "_DEG_down.txt", sep = ""))


###########################################################################
## plot volcano plot
# diffData = fread(file = paste0(outPrefix, "_DEG_all.txt", collapse = ""), sep = "\t", header = T, stringsAsFactors = F, data.table = F)

markGenes <- c()

# tmpDf <- filter(diffData, geneName %in% markGenes)

plotTitle <- paste("Volcano plot:", inter_gt21_BA, sep = " ")
p2 <- volcanoPlot(df = diffData, 
                  namePrefix = outPrefix, 
                  title = plotTitle, 
                  fdr_col = "padj", 
                  lfc_col = "shrinkLog2FC",
                  fdr_Cut = FDR_cut, lfc_cut = lfc_cut, 
                  geneOfInterest = markGenes,
                  ylimit = 10)

# png(filename = paste(outPrefix, "_volcano.png", sep = ""), width = 3000, height = 3000, res = 230)
pdf(file = paste(outPrefix, ".volcano.pdf", sep = ""), width = 10, height = 10)
p2
dev.off()


###########################################################################











library(DESeq2)
library(dplyr)
library(tibble)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tximport)


## This script:
## read the stringtie output using tximport 
## perform differential gene expression analysis using DESeq2


rm(list = ls())
source("E:/Chris_UM/Codes/RNA_Seq/DESeq2_functions.R")

path = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/CAuris_diff"

path_stringtie = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/stringTie/stringTie_CAur"

file_geneInfo = "E:/Chris_UM/Database/Candida_auris_B8441/Candida_auris_B8441.info.tab"
# "E:/Chris_UM/Database/Candida_auris_B8441/Candida_auris_B8441.info.tab"
# "E:/Chris_UM/Database/C_albicans/SC5314_A21/C_albicans_SC5314_A21.info.tab"

file_sampleInfo = "sampleInfo.txt"


outFilePrefix = "CAur_HSP90_depletion_effect"

FDR_cut <- 0.05
lfc_cut = 0.585
up_cut = lfc_cut
down_cut = lfc_cut * -1


path = paste(path, outFilePrefix, sep = "/")

if(!dir.exists(path)){
  dir.create(path = path)
}

setwd(path)



design = ~ genotype + treatment + genotype:treatment

###########################################################################
## set levels for experiment design information

exptInfo = read.table(file = file_sampleInfo, header = T, sep = "\t", stringsAsFactors = F)

## set the reference levels
exptInfo$genotype = factor(exptInfo$genotype, levels = c("WT", "tet90"))
exptInfo$treatment = factor(exptInfo$treatment, levels = c("YPD", "YPD_DOX"))
exptInfo$condition = factor(exptInfo$condition, levels = c("WT_YPD", "WT_YPD_DOX", "tet90_YPD", "tet90_YPD_DOX"))

rownames(exptInfo) = exptInfo$sampleId


###########################################################################
## Add gene information

geneSym = data.table::fread(file = file_geneInfo, sep = "\t", header = T, stringsAsFactors = F) %>%
  distinct(geneId, .keep_all = T)

###########################################################################
## import the counts data using tximport and run DESeq2

filesStringtie = paste(path_stringtie, "/stringTie_", exptInfo$sampleId, "/t_data.ctab", sep = "")
names(filesStringtie) = exptInfo$sampleId

tmp = data.table::fread(file = filesStringtie[1], sep = "\t", header = T, stringsAsFactors = F)
tx2gene = tmp[, c("t_name", "gene_id")]

txi = tximport(files = filesStringtie, type = "stringtie", tx2gene = tx2gene)

ddsTxi = DESeqDataSetFromTximport(txi = txi, colData = exptInfo, design = design)

# assay(ddsTxi)
# colData(ddsTxi)
# rowData(ddsTxi)


## Run DESeq2
dds = DESeq(ddsTxi)


## raw counts
rawCounts = tibble::rownames_to_column(as.data.frame(counts(dds, normalized = FALSE)), var = "geneId")

fwrite(x = rawCounts, file = paste(outFilePrefix, "_rawCounts.tab", sep = ""),
       sep = "\t", row.names = F, col.names = T, quote = F)


## FPKM
fpkmCounts = tibble::rownames_to_column(as.data.frame(fpkm(dds)), var = "geneId")

fwrite(x = fpkmCounts, file = paste(outFilePrefix, "_FPKM.tab", sep = ""),
       sep = "\t", row.names = F, col.names = T, quote = F)


## normalized counts matrix
normCounts = tibble::rownames_to_column(as.data.frame(counts(dds, normalized = TRUE)), var = "geneId")

write.table(x = normCounts, file = paste0(c(outFilePrefix,"_normCounts.tab"), collapse = ""), row.names = F, col.names = T, sep = "\t", quote = F)


## r-log normalized counts
rld = rlog(dds)
rldCount = rownames_to_column(as.data.frame(assay(rld)), var = "geneId")

fwrite(x = rldCount, file = paste(outFilePrefix, "_rlogCounts.tab", sep = ""),
       sep = "\t", row.names = F, col.names = T, quote = F)


###########################################################################
## PCA analysis
##NOTE: Typically, we recommend users to run samples from all groups together, and then use the contrast argument of the results function to extract comparisons of interest after fitting the model using DESeq. The model fit by DESeq estimates a single dispersion parameter for each gene, which defines how far we expect the observed count for a sample will be from the mean value from the model given its size factor and its condition group. However, for some datasets, exploratory data analysis (EDA) plots could reveal that one or more groups has much higher within-group variability than the others. This is case where, by comparing groups A and B separately - subsetting a DESeqDataSet to only samples from those two groups and then running DESeq on this subset - will be more sensitive than a model including all samples together.

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


png(filename = paste(outFilePrefix, "_PCA.png", sep = ""), width = 2500, height = 2000, res = 250)
p1
dev.off()


###########################################################################
## extract results for interaction terms

resultsNames(dds)

## the condition effect B_vs_A for genotype I (the main effect)
res_BA_gt1 = results(dds, name = "treatment_YPD_DOX_vs_YPD")

## another way to extract condition effect for genotype I using contrast argument
# res_BA_gt1 = results(dds, contrast = c("treatment", "YPD_DOX", "YPD"), lfcThreshold = lfc_cut)

summary(res_BA_gt1)
mcols(res_BA_gt1)
plotMA(res_BA_gt1, ylim=c(-4,4), main = "MA plot with unshrunken LFC")


## the condition effect B_vs_A for genotype II
## remember that the contrast has to be a list with first element as a vector of genotype I effect and interaction term
res_BA_gt2 = results(dds,
                     contrast = list(c("treatment_YPD_DOX_vs_YPD", "genotypetet90.treatmentYPD_DOX")))

summary(res_BA_gt2)
mcols(res_BA_gt2)
plotMA(res_BA_gt2, ylim=c(-4,4), main = "MA plot with unshrunken LFC")


## the interaction term for condition effect B_vs_A in genotype II vs genotype I.
## the interaction term, answering: is the condition effect *different* across genotypes
inter_gt21_BA = "genotypetet90.treatmentYPD_DOX"

resInter_gt21_BA = results(dds, name = inter_gt21_BA)

summary(resInter_gt21_BA)
mcols(resInter_gt21_BA)


resShrink_inter_gt21_BA = lfcShrink(dds, coef = inter_gt21_BA,
                                    res = resInter_gt21_BA,
                                    type="apeglm")

summary(resShrink_inter_gt21_BA)
mcols(resShrink_inter_gt21_BA)




###########################################################################
## MA plot

op = par(mfrow = c(2, 1))

plotMA(resInter_gt21_BA, ylim=c(-4,4), main = "MA plot with unshrunken LFC")
plotMA(resShrink_inter_gt21_BA, ylim=c(-4,4), main = "MA plot with apeglm shrunken LFC")

## MA plot with unshrunken LFC
geneplotter::plotMA(object = data.frame(
  m = resInter_gt21_BA$baseMean,
  a = resInter_gt21_BA$log2FoldChange,
  c = ifelse(abs(resInter_gt21_BA$log2FoldChange) >= lfc_cut & resInter_gt21_BA$padj <= FDR_cut, TRUE, FALSE)),
  ylim=c(-4,4), 
  main = "MA plot with unshrunken LFC"
)

## MA plot with apeglm shrunken LFC
geneplotter::plotMA(object = data.frame(
  m = resShrink_inter_gt21_BA$baseMean,
  a = resShrink_inter_gt21_BA$log2FoldChange,
  c = ifelse(abs(resShrink_inter_gt21_BA$log2FoldChange) >= lfc_cut & resShrink_inter_gt21_BA$padj <= FDR_cut, TRUE, FALSE)),
  ylim=c(-4,4), 
  main = "MA plot with apeglm shrunken LFC"
)

par(op)




png(filename = paste0(outFilePrefix, "_MA.png", collapse = ""), width = 2000, height = 3000, res = 250)
op = par(mfrow = c(2, 1))

## MA plot with unshrunken LFC
geneplotter::plotMA(object = data.frame(
  m = resInter_gt21_BA$baseMean,
  a = resInter_gt21_BA$log2FoldChange,
  c = ifelse(abs(resInter_gt21_BA$log2FoldChange) >= lfc_cut & resInter_gt21_BA$padj <= FDR_cut, TRUE, FALSE)),
  ylim=c(-4,4), 
  main = "MA plot with unshrunken LFC"
)

## MA plot with apeglm shrunken LFC
geneplotter::plotMA(object = data.frame(
  m = resShrink_inter_gt21_BA$baseMean,
  a = resShrink_inter_gt21_BA$log2FoldChange,
  c = ifelse(abs(resShrink_inter_gt21_BA$log2FoldChange) >= lfc_cut & resShrink_inter_gt21_BA$padj <= FDR_cut, TRUE, FALSE)),
  ylim=c(-4,4), 
  main = "MA plot with apeglm shrunken LFC"
)

par(op)
dev.off()




###########################################################################
## generate the results DF

## store un-shrunk LFC data
resDf = rownames_to_column(as.data.frame(resInter_gt21_BA), var = "geneId")
head(resDf)

## store shrunk LFC data
resShrinkDf = rownames_to_column(as.data.frame(resShrink_inter_gt21_BA), var = "geneId") 
head(resShrinkDf)


resultTable = resShrinkDf
head(resultTable)

## function to get mean count formula as list
formula_list = function(x){
  ids = samples = sapply(x$sampleId, FUN = as.name, USE.NAMES = F, simplify = T)
  len = nrow(x)
  name = as.name(paste(x$condition[1], "_meanCount", sep = ""))
  
  form = quos(!! name := sum(!!! ids) / !! len)

  return(form)
}


## generate formula list
af = exptInfo %>% 
  # dplyr::mutate_if(.predicate = is.factor, .funs = as.character) %>% 
  dplyr::group_by(condition) %>% 
  dplyr::do(form = formula_list(.)) %>% 
  dplyr::mutate(name = paste(condition, "_meanCount", sep = ""))

avgFormula = unlist(af$form)

## apply formula on each row to count mean count
countsDf = dplyr::rowwise(normCounts) %>%
  dplyr::mutate(!!! avgFormula)

write.table(head(countsDf), file = "clipboard", sep = "\t", row.names = F, quote = F)

selCols = colnames(geneSym)[! colnames(geneSym) %in% c("geneId", "description")]


diffData = resultTable %>% 
  left_join(y = countsDf, by = "geneId") %>% 
  left_join(y = geneSym, by = c("geneId" = "geneId")) %>%
  dplyr::select(geneId, !!!selCols, ends_with("_meanCount"), log2FoldChange, pvalue, padj) %>%
  mutate(
    diff_l2fc = case_when(
      padj < FDR_cut & log2FoldChange >= up_cut ~ "up",
      padj < FDR_cut & log2FoldChange <= down_cut ~ "down",
      TRUE ~ "noDEG"
    )
  )

head(diffData)

significant = filter(diffData, padj < FDR_cut, log2FoldChange >= up_cut | log2FoldChange <= down_cut)

significant_up = filter(diffData, padj < FDR_cut, log2FoldChange >= up_cut)
significant_down = filter(diffData, padj < FDR_cut, log2FoldChange <= down_cut)


write.table(x = resultTable, file = paste0(outFilePrefix, "_DESeq2.txt", collapse = ""),
            col.names = T, row.names = F, sep = "\t", quote = F)
write.table(x = diffData, file = paste0(outFilePrefix, "_DEG_all.txt", collapse = ""),
            col.names = T, row.names = F, sep = "\t", quote = F)
write.table(x = significant_up, file = paste0(outFilePrefix, "_DEG_up.txt", collapse = ""),
            col.names = T, row.names = F, sep = "\t", quote = F)
write.table(x = significant_down, file = paste0(outFilePrefix, "_DEG_down.txt", collapse = ""),
            col.names = T, row.names = F, sep = "\t", quote = F)


###########################################################################
## plot volcano plot
# diffData = fread(file = paste0(outFilePrefix, "_DEG_all.txt", collapse = ""), sep = "\t", header = T, stringsAsFactors = F, data.table = F)

markGenes = c()

# tmpDf = filter(diffData, geneName %in% markGenes)

plotTitle = paste("Volcano plot:", inter_gt21_BA, sep = " ")
p2 = volcanoPlot(df = diffData, 
                 namePrefix = outFilePrefix, 
                 title = plotTitle, 
                 fdr_col = "padj", 
                 lfc_col = "log2FoldChange",
                 fdr_Cut = FDR_cut, lfc_cut = lfc_cut, 
                 geneOfInterest = markGenes,
                 ylimit = 200)

png(filename = paste0(outFilePrefix, "_volcano.png", collapse = ""), width = 3000, height = 3000, res = 230)
p2
dev.off()


###########################################################################











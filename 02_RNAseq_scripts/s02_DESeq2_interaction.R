suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(org.AFumigatus.Af293.eg.db))
suppressPackageStartupMessages(library(here))

## This script:
## read the stringtie output using tximport or as raw counts
## perform differential gene expression analysis with interaction terms using DESeq2

rm(list = ls())
source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

analysisName <- "interaction_znfA_del_casp"

outDir <- here::here("analysis", "02_polII_diff", analysisName)
file_sampleInfo <- here::here("data", "reference_data", "polII_diff_info.txt")

outPrefix <- paste(outDir, analysisName, sep = "/")


orgDb <- org.AFumigatus.Af293.eg.db
key_orgDb <- "GID"

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc
###########################################################################
## set levels for experiment design information

exptInfo <- read.table(file = file_sampleInfo, header = T, sep = "\t", stringsAsFactors = F)

## set the reference levels
exptInfo$genotype <- forcats::as_factor(exptInfo$genotype)
exptInfo$treatment <- forcats::as_factor(exptInfo$treatment)
exptInfo$condition <- forcats::as_factor(exptInfo$condition)

rownames(exptInfo) <- exptInfo$sampleId

design <- ~ genotype + treatment + genotype:treatment

file.copy(from = file_sampleInfo, to = outDir)

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
file_rawCounts <- here::here("data", "polII_data", "polII_raw_counts.tab")

countsDf <- suppressMessages(readr::read_tsv(file = file_rawCounts, col_names = T)) %>%
  as.data.frame()
rownames(countsDf) <- countsDf$geneId
countsDf$geneId <- NULL

## select only those sample rows which are part of current comparison
countsDf <- countsDf[, rownames(exptInfo)]

## run DESeq2 and extract the processed data
ddsCount <- DESeqDataSetFromMatrix(countData = countsDf, colData = exptInfo, design = design)

## Run DESeq2
dds <- DESeq(ddsCount)

###########################################################################
## Add gene symbol for each Ensembl ID
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(x = orgDb, keytype = key_orgDb),
                                  columns = c("GENE_NAME", "DESCRIPTION"),
                                  keytype = key_orgDb) %>% 
  dplyr::rename(geneId = !!key_orgDb)



## Make sure to merge the geneIds which are duplicated using group_by + summarize_all and remove the NA strings in data frame by replacing it with real <NA>
if(any(duplicated(geneInfo$geneId))){
  geneInfo <- dplyr::group_by(geneInfo, geneId) %>%
    summarize_all(.funs = ~ paste0(unique(.), collapse = ","))
  
  ## OR just take first row in case of duplicate rows 
  # geneInfo <- dplyr::distinct(geneInfo, geneId, .keep_all = T)
}

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


pt_pca <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
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


png(filename = paste(outPrefix, ".PCA.png", sep = ""), width = 3000, height = 3000, res = 300)
pt_pca
dev.off()



sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)


pt_dist <- ComplexHeatmap::Heatmap(
  matrix = sampleDistMatrix,
  col = colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255),
  column_title = "Distance matrix of normalized read counts",
  column_title_gp = gpar(fontface = "bold", fontsize = 14),
  heatmap_legend_param = list(title = "Distance", title_gp = gpar(fontsize = 12),
                              title_position = "topcenter")
)


png(filename = paste(outPrefix, ".distance_heatmap.png", sep = ""),
    width = 3000, height = 3000, res = 300)
draw(pt_dist,
     padding = unit(rep(0.5, 4), "cm")
)
dev.off()


###########################################################################
## extract results for interaction terms

resultsNames(dds)

########################
## explained below using following design
## genotype : I | I | I | I | I | I | II | II | II | II | II | II
## condition: A | A | A | B | B | B |  A |  A |  A |  B |  B |  B
## ~ genotype + condition + genotype:condition
########################

########################
## the condition effect B_vs_A for genotype I (the main effect)
## contrast=c("condition","B","A")
res_BA_gt1 <- results(dds, name = "treatment_casp_vs_ctrl")

summary(res_BA_gt1)
mcols(res_BA_gt1)

resLfcShrink_BA_gt1 <- lfcShrink(dds, coef = "treatment_casp_vs_ctrl",
                                 res = res_BA_gt1,
                                 type="apeglm")

op <- par(mfrow = c(2, 1))
plotMA(res_BA_gt1, ylim=c(-4,4), main = "MA plot with unshrunken LFC")
plotMA(resLfcShrink_BA_gt1, ylim=c(-4,4), main = "MA plot with shrunken LFC")
par(op)

readr::write_tsv(x = rownames_to_column(as.data.frame(res_BA_gt1), var = "geneId"),
                 path = paste(outPrefix, ".casp_vs_ctrl.WT.DESeq2.tab", sep = ""))


########################
## the condition effect B_vs_A for genotype II:
## this is, by definition, the main effect *plus* the interaction term
## (the extra condition effect in genotype II compared to genotype I).
## results(dds, list( c("condition_B_vs_A","genotypeII.conditionB") ))
##
## remember that the contrast has to be a list with first element as a vector of 
## genotype I effect and interaction term

res_BA_gt2 <- results(
  dds,
  contrast = list(c("treatment_casp_vs_ctrl", "genotypedel_znfA.treatmentcasp"))
)

summary(res_BA_gt2)
mcols(res_BA_gt2)

## contrast cannot be used with apeglm. so use ashr for shrinking LFC
resLfcShrin_BA_gt2 <- lfcShrink(
  dds, res = res_BA_gt2, type = "ashr",
  contrast = list(c("treatment_casp_vs_ctrl", "genotypedel_znfA.treatmentcasp"))
)

op <- par(mfrow = c(2, 1))
plotMA(res_BA_gt2, ylim=c(-4,4), main = "MA plot with unshrunken LFC")
plotMA(resLfcShrin_BA_gt2, ylim=c(-4,4), main = "MA plot with shrunken LFC")
par(op)

readr::write_tsv(x = rownames_to_column(as.data.frame(res_BA_gt2), var = "geneId"),
                 path = paste(outPrefix, ".casp_vs_ctrl.del_znfA.DESeq2.tab", sep = ""))

########################
## the interaction term for condition effect B_vs_A in genotype II vs genotype I.
## the interaction term, answering: is the condition effect *different* across genotypes
inter_gt21_BA <- "genotypedel_znfA.treatmentcasp"

resInter_gt21_BA <- results(dds, name = inter_gt21_BA)

summary(resInter_gt21_BA)
mcols(resInter_gt21_BA)

resultSummary <- paste(capture.output(summary(resInter_gt21_BA))[1:8], collapse = "\n")

resShrink_inter_gt21_BA <- lfcShrink(dds, coef = inter_gt21_BA,
                                     res = resInter_gt21_BA,
                                     type="apeglm")

summary(resShrink_inter_gt21_BA)
mcols(resShrink_inter_gt21_BA)

###########################################################################

fcDensity <- hist(x = pmin(pmax(resInter_gt21_BA$log2FoldChange, -5), 5),
                  breaks = 50,
                  main = "log2(fold-change) frequency distribution")

hist(x = resInter_gt21_BA$pvalue, breaks = 100, main = "p-value distribution")
hist(x = resInter_gt21_BA$pvalue, breaks = c(seq(0, 0.1, length.out = 20), seq(0.11,1, by = 0.01)),
     main = "p-value distribution")

hist(x = resInter_gt21_BA$padj, breaks = 100, main = "q-value distribution")

pt_lfcFreq <- ggplot() +
  geom_histogram(
    data = as.data.frame(resInter_gt21_BA),
    mapping = aes(x = pmin(pmax(log2FoldChange, -5), 5), y = ..density.. , fill = "unshrunken"),
    bins = 100, alpha = 0.5) +
  geom_histogram(
    data = as.data.frame(resShrink_inter_gt21_BA),
    mapping = aes(x = pmin(pmax(log2FoldChange, -5), 5), y = ..density.. , fill = "shrunken"),
    bins = 100, alpha = 0.5) +
  annotate(geom = "text", x = -4.5, y = Inf, label = resultSummary, vjust = 1, hjust = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  scale_x_continuous(limits = c(-5, 5), expand = expand_scale(mult = 0.01)) +
  scale_y_continuous(expand = expand_scale(mult = 0.01)) +
  scale_fill_manual(
    name = NULL,
    values = c("unshrunken" = "#999999", "shrunken" = "red"),
    breaks = c("unshrunken", "shrunken"),
    labels = c("unshrunken LFC", "shrunken LFC")
  ) +
  labs(title = paste("log2(fold change) frequency distribution:", analysisName),
       x = "log2(fold change)", y = "Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(face = "bold", size = 15),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.2,"cm")
  )


pt_pqDensity <- ggplot(data = as.data.frame(resInter_gt21_BA)) +
  geom_histogram(mapping = aes(x = pvalue, y = ..density.. , fill = "pvalue"),
                 bins = 100, alpha = 0.5) +
  geom_histogram(mapping = aes(x = padj, y = ..density.. , fill = "padj"),
                 bins = 100, alpha = 0.5) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
  annotate(geom = "text", x = 0.3, y = Inf, label = resultSummary, vjust = 1) +
  scale_fill_manual(
    name = NULL,
    values = c("pvalue" = "#999999", "padj" = "#E69F00"),
    breaks = c("pvalue", "padj"),
    labels = c("p-value", "q-value")
  ) +
  scale_x_continuous(expand = expand_scale(mult = 0.01)) +
  scale_y_continuous(expand = expand_scale(mult = 0.01)) +
  labs(title = paste("p-value and q-value density distribution:", analysisName),
       x = "p-value or q-value", y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.2,"cm"),
        # legend.direction = "horizontal",
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.title = element_text(face = "bold", size = 15),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1)
  )

###########################################################################
## MA plot
op <- par(mfrow = c(2, 1))

plotMA(resInter_gt21_BA, ylim=c(-4,4), main = "MA plot with unshrunken LFC")
plotMA(resShrink_inter_gt21_BA, ylim=c(-4,4), main = "MA plot with apeglm shrunken LFC")

par(op)


###########################################################################
## generate the results DF

## store un-shrunk LFC data
resDf <- rownames_to_column(as.data.frame(resInter_gt21_BA), var = "geneId")
glimpse(resDf)

## store shrunk LFC data
resShrinkDf <- rownames_to_column(as.data.frame(resShrink_inter_gt21_BA), var = "geneId") 
glimpse(resShrinkDf)

resultTable <- dplyr::left_join(
  x = resDf,
  y = dplyr::select(resShrinkDf, geneId, shrinkLog2FC = log2FoldChange, sLfcSE = lfcSE),
  by = c("geneId")) %>%
  dplyr::select(geneId, baseMean, log2FoldChange, lfcSE, shrinkLog2FC, sLfcSE, everything())

# resultTable <- resShrinkDf
glimpse(resultTable)



## function to get mean count formula as list
formula_list <- function(x){
  ids = samples = sapply(x$sampleId, FUN = as.name, USE.NAMES = F, simplify = T)
  len = nrow(x)
  name = as.name(paste(x$condition[1], ".meanCount", sep = ""))
  
  form = quos(!! name := sum(!!! ids) / !! len)
  
  return(form)
}


## generate formula list
af <- exptInfo %>% 
  # dplyr::mutate_if(.predicate = is.factor, .funs = as.character) %>% 
  dplyr::group_by(condition) %>% 
  dplyr::do(form = formula_list(.)) %>% 
  dplyr::mutate(name = paste(condition, ".meanCount", sep = ""))

avgFormula <- unlist(af$form)

## apply formula on each row to count mean count
countsDf <- dplyr::rowwise(normCounts) %>%
  dplyr::mutate(!!! avgFormula)


selCols <- colnames(geneInfo)[colnames(geneInfo) != "geneId"]

diffData <- resultTable %>% 
  left_join(y = countsDf, by = "geneId") %>% 
  left_join(y = geneInfo, by = c("geneId" = "geneId")) %>%
  mutate(
    diff_l2fc = case_when(
      padj < cutoff_fdr & log2FoldChange >= cutoff_up ~ "up",
      padj < cutoff_fdr & log2FoldChange <= cutoff_down ~ "down",
      TRUE ~ "noDEG"
    )
  ) %>% 
  dplyr::select(geneId, ends_with(".meanCount"), log2FoldChange, shrinkLog2FC,
                pvalue, padj, diff_l2fc, !!!selCols)

glimpse(diffData)

significant <- filter(diffData, padj < cutoff_fdr, log2FoldChange >= cutoff_up | log2FoldChange <= cutoff_down)

significant_up <- filter(diffData, padj < cutoff_fdr, log2FoldChange >= cutoff_up)
significant_down <- filter(diffData, padj < cutoff_fdr, log2FoldChange <= cutoff_down)

readr::write_tsv(x = resultTable, path = paste(outPrefix, ".DESeq2.tab", sep = ""))
readr::write_tsv(x = diffData, path = paste(outPrefix, ".DEG_all.txt", sep = ""))

# readr::write_tsv(x = significant_up, path = paste(outPrefix, "_DEG_up.txt", sep = ""))
# readr::write_tsv(x = significant_down, path = paste(outPrefix, "_DEG_down.txt", sep = ""))


###########################################################################
## plot volcano plot
# diffData = fread(file = paste0(outPrefix, "_DEG_all.txt", collapse = ""), sep = "\t", header = T, stringsAsFactors = F, data.table = F)

markGenes <- c()

plotTitle <- paste("Volcano plot:", inter_gt21_BA, sep = " ")
pt_volc <- volcano_plot(df = diffData,
                   title = plotTitle,
                   fdr_col = "padj",
                   lfc_col = "log2FoldChange",
                   fdr_cut = cutoff_fdr, lfc_cut = cutoff_lfc, 
                   markGenes = markGenes,
                   ylimit = 4, xlimit = c(-4, 4))

png(filename = paste(outPrefix, ".volcano.png", sep = ""), width = 2500, height = 3000, res = 280)
plot(pt_volc$plot)
dev.off()


###########################################################################



# plot all data in single PDF file

pdf(file = paste(outPrefix, ".summary_plots.pdf", sep = ""), width = 10, height = 10, onefile = TRUE)

## PCA
plot(pt_pca)

## distance heatmap
draw(pt_dist,
     padding = unit(rep(0.5, 4), "cm")
)

## MA plots
op <- par(mfrow = c(2, 1))
plotMA(
  object = resInter_gt21_BA,
  ylim=c(-4,4),
  main = paste("MA plot with unshrunken LFC:", analysisName)
)

plotMA(
  object = resShrink_inter_gt21_BA,
  ylim=c(-4,4),
  main = paste("MA plot with shrunken LFC:", analysisName)
)

par(op)

## volcano plot
plot(pt_volc$plot)

## p-value distribution plots
plot(pt_pqDensity)
plot(pt_lfcFreq)

dev.off()

###########################################################################











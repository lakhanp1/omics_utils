suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(require(openxlsx))


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

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################
## configuration and cutoffs

diffDataPath <- here::here("analysis", "07_polII_diff")
file_sampleInfo <- here::here("data", "reference_data", "polII_sample_info.txt")

useAllGroupsSamples <- FALSE

cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- cutoff_lfc * -1

orgDb <- org.Anidulans.FGSCA4.eg.db
col_geneId <- "GID"
otherCols <- c("GENE_NAME", "DESCRIPTION")

###########################################################################
#############################################
## Run DESeq2 pipeline using Rscript       ##
## command line arguments                  ##
#############################################

parser <- ArgumentParser(
  description = "This script automates the pairwise differential gene expression analysis by DESeq2. It reads a information from tabular config file for a given comparison."
)

parser$add_argument(
  "-c", "--config", action="store",
  dest = "config", type = "character", nargs = 1, required = TRUE,
  # default = here::here("data", "reference_data", "polII_DESeq2_info.txt"),
  help = "DEG configuration TAB delimited file with columns: comparison, type, group1, group2, design, samples"
)

parser$add_argument(
  "-d", "--deg", action="store",
  dest = "deg", type = "character", nargs = 1, required = TRUE,
  help = "DEG comparison ID. This value should be present in the column comparison of config file"
)

# parser$print_help()

# file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
# analysisName <- "AN10021_sCopy_OE_vs_MH11036"

args <- parser$parse_args()

file_RNAseq_info <- args$config
analysisName <- args$deg


rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>%
  dplyr::filter(comparison == analysisName, type == "pairwise")

if(nrow(rnaseqInfo) != 1){
  stop(analysisName, " RNAseq data does not exist in config file")
}

## the denominator or WT in log2(fold_change) should be second
compare <- c(rnaseqInfo$group1, rnaseqInfo$group2)
col_compare <- rnaseqInfo$design
samples <- unlist(stringr::str_split(string = rnaseqInfo$samples, pattern = ";"))

# #############################################
# ## User input of conditions to be compared ##
# #############################################
# analysisName <- "SCX1_KO_ctrl_vs_WT_ctrl"
# 
# ## the denominator or WT in log2(fold_change) should be second
# compare <- c("SCX1_KO_ctrl", "SCX1_WT_ctrl")
# col_compare <- "conditionSC"


##########################################################################

outDir <- paste(diffDataPath, analysisName, sep = "/")
outPrefix <- paste(outDir, analysisName, sep = "/")

if(!dir.exists(outDir)){
  dir.create(path = outDir, recursive = TRUE)
}

exptInfo <- read.table(file = file_sampleInfo, header = T, sep = "\t", stringsAsFactors = F)

## get subset of experiment information df
if(isFALSE(useAllGroupsSamples)){
  # exptInfo <- dplyr::filter(exptInfo, condition %in% compare)
  exptInfo <- dplyr::filter(exptInfo, sampleId %in% samples)
  
  if(nrow(exptInfo) == 0){
    stop("No samples in the experiment information")
  }
}

if(isFALSE(setequal(compare, exptInfo[[col_compare]]))){
  stop("Mismatch in the ", col_compare, "s to compare and values from ", col_compare, " column of exptInfo")
}


# ## ensure that reference level is same as compare[2] in the factor
exptInfo <- exptInfo %>% dplyr::mutate(
  !!col_compare := forcats::fct_relevel(.f = !!sym(col_compare), compare[2], compare[1])
)

rownames(exptInfo) <- exptInfo$sampleId

# design <- ~ condition
design <- as.formula(paste("~", col_compare))

###########################################################################
## import counts data: either by tximport or as raw count matrix

## import the counts data using tximport and run DESeq2
readLength <- as.numeric(readr::read_file(file = here::here("data", "read_length.config")))
path_stringtie <- here::here("data", "stringTie")
filesStringtie <- paste(path_stringtie, "/stringTie_", exptInfo$stringtieId, "/t_data.ctab", sep = "")
names(filesStringtie) <- exptInfo$sampleId

tmp <- data.table::fread(file = filesStringtie[1], sep = "\t", header = T, stringsAsFactors = F)
tx2gene <- tmp[, c("t_name", "gene_id")]

txi <- tximport(files = filesStringtie, type = "stringtie",
                tx2gene = tx2gene, readLength = readLength)

ddsTxi <- DESeqDataSetFromTximport(txi = txi, colData = exptInfo, design = design)
# assay(ddsTxi)
# colData(ddsTxi)
# rowData(ddsTxi)

## Run DESeq2
dds <- DESeq(ddsTxi)


# ## import raw counts data and run DESeq2
# file_rawCounts <- here::here("data", "polII_data", "polII_raw_counts.tab")
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
# # dplyr::glimpse(countsDf)
# 
# ## run DESeq2 and extract the processed data
# ddsCount <- DESeqDataSetFromMatrix(countData = countsDf, colData = exptInfo, design = design)
# 
# ## Run DESeq2
# dds <- DESeq(ddsCount)

###########################################################################
## Add gene symbol for each Ensembl ID
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(x = orgDb, keytype = col_geneId),
                                  columns = otherCols,
                                  keytype = col_geneId) %>%
  dplyr::rename(geneId = !!col_geneId)



## Make sure to merge the geneIds which are duplicated using group_by + summarize_all and remove the NA strings in data frame by replacing it with real <NA>
if(any(duplicated(geneInfo$geneId))){
  geneInfo <- dplyr::group_by(geneInfo, geneId) %>%
    summarize_all(.funs = ~ paste0(unique(.), collapse = ","))
  
  ## OR just take first row in case of duplicate rows
  # geneInfo <- dplyr::distinct(geneInfo, geneId, .keep_all = T)
}

# dplyr::glimpse(geneInfo)

###########################################################################

## raw counts
rawCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = FALSE)), var = "geneId")
readr::write_tsv(x = rawCounts, path = paste(outPrefix, ".rawCounts.tab", sep = ""))

# ## FPKM
# fpkmCounts <- tibble::rownames_to_column(as.data.frame(fpkm(dds)), var = "geneId")
# readr::write_tsv(x = fpkmCounts, path = paste0(c(outPrefix,".FPKM.tab"), collapse = ""))

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


plotPCA(rld, intgroup=c(col_compare))

pcaData <- plotPCA(rld, intgroup=c(col_compare), returnData = TRUE) %>% 
  tibble::rownames_to_column(var = "sampleId") %>% 
  dplyr::left_join(y = exptInfo, by = c("sampleId", col_compare))

percentVar <- sprintf("%.2f", 100 * attr(pcaData, "percentVar"))

pltTitle <- paste("Principal Component Analysis:", compare[1], "vs", compare[2])

pcaData <- dplyr::mutate(
  pcaData,
  # !!sym(treatment) := forcats::as_factor(!!sym(treatment)),
  # !!sym(time) := forcats::as_factor(!!sym(time)),
  !!sym(col_compare) := forcats::as_factor(!!sym(col_compare))
)

fillColumn <- col_compare
# shapeColumn <- "time"

if(length(unique(pcaData[[fillColumn]])) <= 9){
  pointCol <- base::structure(
    .Data = RColorBrewer::brewer.pal(n = length(unique(pcaData[[fillColumn]])), name = "Set1"),
    names = levels(pcaData[[fillColumn]])
  )
} else{
  pointCol <- base::structure(
    .Data = rainbow(n = length(unique(pcaData[[fillColumn]]))),
    names = levels(pcaData[[fillColumn]])
  )
}

pt_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.2,"cm"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.title = element_text(face = "bold", size = 15)
  )

pt_pca <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = !!sym(fillColumn)),
             size = 6, stroke = 2) +
  # scale_shape_manual(values = c(1, 15, 17)) +
  # guides(fill=guide_legend(override.aes=list(shape=21))) +
  scale_color_manual(values = pointCol) +
  geom_text_repel(mapping = aes(label = name), size = 4,  
                  point.padding = unit(0.5, "lines")) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle(pltTitle) +
  pt_theme +
  theme(
    legend.position = "bottom",
  )


# png(filename = paste(outPrefix, ".PCA.png", sep = ""), width = 3000, height = 3000, res = 300)
# pt_pca
# dev.off()

## PCA on subset of the data
# subExptInfo <- dplyr::filter(exptInfo, condition %in% compare)
#
# seTemp <- SummarizedExperiment(
#   assays = log2(counts(dds, normalized = TRUE) + 1)[, subExptInfo$sampleId],
#   colData = subExptInfo
# )
#
# plotPCA( DESeqTransform(seTemp))

###########################################################################
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)


pt_dist <- ComplexHeatmap::Heatmap(
  matrix = sampleDistMatrix,
  col = colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255),
  row_names_gp = gpar(fontface = "bold", fontsize = 18),
  column_names_gp = gpar(fontface = "bold", fontsize = 18),
  row_names_max_width = max_text_width(rownames(sampleDistMatrix), gp = gpar(fontsize = 18)),
  column_names_max_height = max_text_width(rownames(sampleDistMatrix), gp = gpar(fontsize = 18)),
  heatmap_legend_param = list(
    title = "Distance", title_gp = gpar(fontsize = 16, fontface = "bold"),
    title_position = "lefttop", direction = "horizontal",
    labels_gp = gpar(fontsize = 16)
  )
)


###########################################################################
## result for the comparison of interest
## "Independent Filtering": For weakly expressed genes, we have no chance of seeing differential expression, because the low read counts suffer from such high Poisson noise that any biological effect is drowned in the uncertainties from the sampling at a low rate. Genes with very low mean count have little or no power, and are best excluded from testing. At first sight, there may seem to be little benefit in filtering out these genes. After all, the test found them to be non-significant anyway. However, these genes have an influence on the multiple testing adjustment, whose performance improves if such genes are removed. By removing the low count genes from the input to the FDR procedure, we can find more genes to be significant among those that we keep, and so improved the power of our test. This approach is known as independent filtering. REF: https://www.bioconductor.org/help/workflows/rnaseqGene/#independent-filtering

resultsNames(dds)
# coefName <- "condition_5A9_Control_vs_CEA17_Control"
coefName <- paste(col_compare, "_", compare[1], "_vs_", compare[2], sep = "")

if(!any(coefName %in% resultsNames(dds))){
  stop("Wrong coefficient name")
}

res <- results(dds, name = coefName, alpha = 0.05)
# summary(res)

contrast <- stringr::str_replace(
  string = mcols(res)["log2FoldChange", ]$description,
  pattern = ".*: ", replacement = "") %>% 
  stringr::str_replace_all(c(" " = "_"))

## get shrunken log2 fold changes
resShrink <- lfcShrink(dds, coef = coefName,
                       res=res, type="apeglm")

resultSummary <- paste(capture.output(summary(resShrink))[1:8], collapse = "\n")

readr::write_lines(
  x = c(paste("DESeq2 analysis:", analysisName), resultSummary),
  path = paste(outPrefix, ".DESeq2_summary.txt", sep = ""))

# mcols(resShrink, use.names=TRUE)


fcDensity <- hist(x = pmin(pmax(res$log2FoldChange, -5), 5),
                  breaks = 50,
                  main = "log2(fold-change) frequency distribution")

# hist(x = res$pvalue, breaks = 100, main = "p-value distribution")
# hist(x = res$pvalue, breaks = c(seq(0, 0.1, length.out = 20), seq(0.11,1, by = 0.01)),
#      main = "p-value distribution")
# 
# hist(x = res$padj, breaks = 100, main = "q-value distribution")


pt_theme <- pt_theme +
  theme(legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1))

## Independent filtering plot
pt_indFil <- ggplot(data = metadata(res)$filterNumRej,
                    mapping = aes(x = theta, y = numRej)) +
  geom_point(shape = 1, size = 3, stroke = 1.5) +
  geom_line(size = 1) +
  geom_line(data = metadata(res)$lo.fit %>% as_tibble(),
            mapping = aes(x = x, y = y),
            color = "red", size = 1) +
  geom_vline(xintercept = metadata(res)$filterTheta,
             color = "blue", linetype = "dashed", size = 1) +
  scale_x_continuous(
    labels = scales::percent_format(), breaks = seq(0, 1, length.out = 6)) +
  labs(
    title = "Independent filtering cutoff by DESeq2",
    subtitle = paste(
      "alpha = ", metadata(res)$alpha, " || filter threshold quantile = ", 
      scales::percent(metadata(res)$filterTheta, accuracy = 0.01),
      "; || filter threshold read count = ", round(metadata(res)$filterThreshold, 2), sep = ""),
    x = "quantiles of filter",
    y = "genes with adjusted p-value < 0.05"
  ) +
  pt_theme


## log2FoldChange density plot
pt_lfcFreq <- ggplot() +
  geom_histogram(
    data = as.data.frame(res),
    mapping = aes(x = pmin(pmax(log2FoldChange, -5), 5), y = ..density.. , fill = "unshrunken"),
    bins = 100, alpha = 0.5) +
  geom_histogram(
    data = as.data.frame(resShrink),
    mapping = aes(x = pmin(pmax(log2FoldChange, -5), 5), y = ..density.. , fill = "shrunken"),
    bins = 100, alpha = 0.5) +
  annotate(geom = "text", x = -4.5, y = Inf, label = resultSummary, vjust = 1, hjust = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +
  scale_x_continuous(limits = c(-5, 5), expand = expansion(mult = 0.01)) +
  scale_y_continuous(expand = expansion(mult = 0.01)) +
  scale_fill_manual(
    name = NULL,
    values = c("unshrunken" = "#999999", "shrunken" = "red"),
    breaks = c("unshrunken", "shrunken"),
    labels = c("unshrunken LFC", "shrunken LFC")
  ) +
  labs(title = paste("log2(fold change) frequency distribution:", compare[1], "vs", compare[2]),
       x = "log2(fold change)", y = "Frequency") +
  pt_theme


## pvalue and qvalue density plot
pt_pqDensity <- ggplot(data = as.data.frame(res)) +
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
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  scale_y_continuous(expand = expansion(mult = 0.01)) +
  labs(title = paste("p-value and q-value density distribution:", compare[1], "vs", compare[2]),
       x = "p-value or q-value", y = "Density") +
  pt_theme


###########################################################################
## MA plot
# png(filename = paste(outPrefix, ".MA.png", sep = ""), width = 2000, height = 3000, res = 250)
# pdf(file = paste(outPrefix, ".MA.pdf", sep = ""), width = 8, height = 10, onefile = F)
op <- par(mfrow = c(2, 1))

plotMA(res, ylim=c(-4,4), main = "MA plot with unshrunken LFC")
plotMA(resShrink, ylim=c(-4,4), main = "MA plot with shrunken log2 fold changes")

par(op)
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
grp1 <- sapply(
  X = rownames(exptInfo[exptInfo[[col_compare]] %in% compare[1], ]),
  FUN = as.name, USE.NAMES = F, simplify = T
)
name1 <- paste(compare[1], "_meanCount", sep = "")
# name1 <- "condition1_meanCount"
grp1Len <- length(grp1)

grp2 <- sapply(
  X = rownames(exptInfo[exptInfo[[col_compare]] %in% compare[2], ]),
  FUN = as.name, USE.NAMES = F, simplify = T
)
name2 <- paste(compare[2], "_meanCount", sep = "")
# name2 <- "condition2_meanCount"
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
      padj < cutoff_fdr & log2FoldChange >= cutoff_up ~ "up",
      padj < cutoff_fdr & log2FoldChange <= cutoff_down ~ "down",
      TRUE ~ "noDEG"
    ),
    # diff_shrink_l2fc = dplyr::case_when(
    #   padj < cutoff_fdr & shrinkLog2FC >= cutoff_up ~ "up",
    #   padj < cutoff_fdr & shrinkLog2FC <= cutoff_down ~ "down",
    #   TRUE ~ "noDEG"
    # )
  )



degData <- diffData %>%
  dplyr::select(
    geneId, contrast, ends_with("_meanCount"), log2FoldChange,
    shrinkLog2FC, pvalue, padj, diff_l2fc, !!!selCols)

# dplyr::glimpse(degData)

significant_up <- filter(degData, padj < cutoff_fdr, log2FoldChange >= cutoff_up)
significant_down <- filter(degData, padj < cutoff_fdr, log2FoldChange <= cutoff_down)

###########################################################################
## store data

readr::write_tsv(x = resDf, path = paste(outPrefix, ".DESeq2.tab", sep = ""))
readr::write_tsv(x = resShrinkDf, path = paste(outPrefix, ".DESeq2_shrunken.tab", sep = ""))
readr::write_tsv(x = degData, path = paste(outPrefix, ".DEG_all.txt", sep = ""))


## write data to excel file
wb <- openxlsx::createWorkbook(creator = "Lakhansing Pardehi Genomics Core")
openxlsx::addWorksheet(wb = wb, sheetName = "DESeq2_DEGs")
openxlsx::writeData(
  wb = wb, sheet = 1, startCol = 2, startRow = 1,
  x = paste("## Differential gene expression analysis by DESeq2 for",
            col_compare,":", compare[1], "/", compare[2],
            "## log2FoldChange cutoff =", cutoff_up, "(up) /", cutoff_down, "(down);",
            "q-value cutoff =", cutoff_fdr)
)
openxlsx::writeData(
  wb = wb, sheet = 1, x = degData,
  startCol = 1, startRow = 2, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
headerStyle <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#e6e6e6")
openxlsx::addStyle(wb = wb, sheet = 1, style = headerStyle, rows = 2, cols = 1:ncol(degData))
openxlsx::setColWidths(wb = wb, sheet = 1, cols = 1, widths = "auto")
openxlsx::setColWidths(wb = wb, sheet = 1, cols = 2, widths = str_length(coefName))
openxlsx::freezePane(wb = wb, sheet = 1, firstActiveRow = 3, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(wb = wb, file = paste(outPrefix, ".DEG_all.xlsx", sep = ""), overwrite = TRUE)

###########################################################################
## plot volcano plot

markGenes <- c()

plotTitle <- paste("Volcano plot:", compare[1], "vs", compare[2], sep = " ")
pt_volc <- volcano_plot(df = diffData,
                        title = plotTitle,
                        fdr_col = "padj",
                        lfc_col = "log2FoldChange",
                        fdr_cut = cutoff_fdr, lfc_cut = cutoff_lfc, 
                        markGenes = markGenes,
                        ylimit = 4, xlimit = c(-4, 4))

# png(filename = paste(outPrefix, ".volcano.png", sep = ""), width = 2500, height = 3000, res = 280)
# plot(pt_volc$plot)
# dev.off()

###########################################################################

pt_pca_volc <- ggpubr::ggarrange(pt_pca, pt_volc$plot, ncol = 2, align = "h")
png(filename = paste(outPrefix, ".PCA_volcano.png", sep = ""), width = 5000, height = 3000, res = 300)
print(pt_pca_volc)
dev.off()

# plot all data in single PDF file

pdf(file = paste(outPrefix, ".summary_plots.pdf", sep = ""), width = 10, height = 10, onefile = TRUE)

## PCA
print(pt_pca)

## distance heatmap
draw(
  pt_dist,
  column_title = paste("Distance matrix of normalized read counts:", analysisName),
  column_title_gp = gpar(fontface = "bold", fontsize = 16),
  padding = unit(rep(0.5, 4), "cm"),
  heatmap_legend_side = "bottom"
)

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

print(pt_indFil)

## volcano plot
print(pt_volc$plot)

## p-value distribution plots
print(pt_lfcFreq)
print(pt_pqDensity)

dev.off()


###########################################################################

# ## heatmap of normalized read counts
# significant <- filter(diffData, padj < cutoff_fdr, log2FoldChange >= cutoff_up | log2FoldChange <= cutoff_down)
# 
# ## rld z-score heatmap
# geneCounts <- significant %>%
#   dplyr::select(geneId, !!!c(grp1, grp2)) %>%
#   dplyr::distinct() %>%
#   column_to_rownames(var = "geneId")
# 
# countMat <- data.matrix(geneCounts)
# 
# countZscoreMat <- chipmine::scale_matrix_rows(x = countMat)
# 
# ## plot main rld score heatmap
# countHeatmap <- Heatmap(
#   countZscoreMat,
#   name = "count_heatmap",
#   col = colorRamp2(breaks = c(min(countZscoreMat), 0, max(countZscoreMat)),
#                    colors = c("green", "black", "red"), space = "LAB"),
#   show_row_names = FALSE,
#   column_names_gp = gpar(fontsize = 14),
#   cluster_columns = FALSE,
#   width = unit(10, "cm"), row_names_max_width = unit(15, "cm"),
#   heatmap_legend_param = list(title = "z-score(normalized read count)", color_bar = "continuous")
# )
# 
# 
# ## fold change heatmap
# foldChangeDf <- dplyr::select(significant, geneId, log2FoldChange) %>%
#   dplyr::distinct() %>%
#   tibble::column_to_rownames(var = "geneId")
# 
# if(all(rownames(geneCounts) != rownames(foldChangeDf))){
#   stop("gene order does not match")
# }
# 
# foldChangeMat <- data.matrix(foldChangeDf)
# 
# 
# fcHeatmap <- Heatmap(matrix = foldChangeMat,
#                      name = "lfc_heatmap",
#                      col = colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"), space = "LAB"),
#                      cluster_rows = TRUE,
#                      clustering_distance_rows = "euclidean",
#                      cluster_columns = TRUE,
#                      show_row_names = FALSE,
#                      column_names_gp = gpar(fontsize = 14),
#                      width = unit(2, "cm"),
#                      heatmap_legend_param = list(title = "\nlog2(fold_change)")
# )
# 
# htList <- countHeatmap + fcHeatmap
# 
# pdf(file = paste(outPrefix, ".normCount_heatmap.pdf", sep = ""), width = 10, height = 12, onefile = TRUE)
# 
# draw(object = htList,
#      main_heatmap = "lfc_heatmap",
#      column_title = paste("Normalized read counts heatmap:", analysisName),
#      row_title = "Genes",
#      row_dend_side = "left",
#      column_title_gp = gpar(fontsize = 14, fontface = "bold")
# )
# 
# dev.off()
# 
# ###########################################################################














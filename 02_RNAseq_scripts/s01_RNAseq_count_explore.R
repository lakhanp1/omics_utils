# suppressPackageStartupMessages(library(ballgown))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(PoiClaClu))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(org.HSapiens.gencodev30.eg.db))

## This script uses tximport to extract the FPKM matrix from the stringTie output
## It also perform PCA analysis using FPKM values
rm(list = ls())

analysisName <- "count_data"

file_sampleInfo <- here::here("data", "reference_data", "sample_info.txt")
readLength <- as.numeric(readr::read_file(file = here::here("data", "reference_data", "read_length.config")))

outDir <- here("analysis", "01_count_data")
outPrefix <- paste(outDir, "/", analysisName, sep = "")
orgDb <- org.HSapiens.gencodev30.eg.db

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


###########################################################################
exptInfo <- suppressMessages(readr::read_tsv(file = file_sampleInfo)) %>% 
  as.data.frame()

rownames(exptInfo) <- exptInfo$sampleId

## set the factor levels.
## control levels should be first
# exptInfo$gt <- forcats::as_factor(exptInfo$gt)
# exptInfo$treatment <- forcats::as_factor(exptInfo$treatment)
exptInfo$condition <- forcats::as_factor(exptInfo$condition)

## Add gene symbol for each Ensembl ID
geneInfo <- AnnotationDbi::select(
  x = orgDb, keytype = "ENSEMBL_VERSION", 
  keys = keys(x = orgDb, keytype = "ENSEMBL_VERSION"),
  columns = c("ENSEMBL", "GENE_NAME", "DESCRIPTION")
) %>% 
  dplyr::rename(geneId = ENSEMBL_VERSION)

design <- ~ condition

###########################################################################
## import counts data: either by tximport or as raw count matrix

## import the counts data using tximport and run DESeq2
path_stringtie <- here::here("data", "stringTie")
filesStringtie <- paste(path_stringtie, "/stringTie_", exptInfo$stringtieId, "/t_data.ctab", sep = "")
names(filesStringtie) <- exptInfo$sampleId

tmp <- data.table::fread(file = filesStringtie[1], sep = "\t", header = T, stringsAsFactors = F)
tx2gene <- tmp[, c("t_name", "gene_id")]

cat("Importing stringTie data with tximport. Read length = ", readLength)

txi <- tximport(files = filesStringtie, type = "stringtie",
                tx2gene = tx2gene, readLength = readLength)

ddsTxi <- DESeqDataSetFromTximport(txi = txi, colData = exptInfo, design = design)
# assay(ddsTxi)
# colData(ddsTxi)
# rowData(ddsTxi)

## Run DESeq2
dds <- DESeq(ddsTxi)

# ## import raw counts data and run DESeq2
# file_rawCounts <- here::here("data", "MatrixCountsPerGeneBySample.Reneto.tab")
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
## raw counts
rawCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = FALSE)), var = "geneId")
readr::write_tsv(x = rawCounts, path = paste0(c(outPrefix,".rawCounts.tab"), collapse = ""))

## FPKM
fpkmCounts <- tibble::rownames_to_column(as.data.frame(fpkm(dds)), var = "geneId");

if(all(fpkmCounts$geneId %in% geneInfo$geneId)){
  fpkmCounts <- dplyr::left_join(x = fpkmCounts, y = geneInfo, by = c("geneId" = "geneId")) %>% 
    dplyr::select(geneId, exptInfo$sampleId, everything()) %>% 
    dplyr::filter(!is.na(geneId))
  
  readr::write_tsv(x = fpkmCounts, path = paste0(c(outPrefix,".FPKM.tab"), collapse = ""))
  
} else{
  warning(
    "Missing geneIds in fpkmCounts from org.db:\n",
    paste(head(fpkmCounts$geneId[which(!fpkmCounts$geneId %in% geneInfo$geneId)], 10), collapse = ", "),
    ", ..."
  )
}


## normalized counts matrix
normCounts <- tibble::rownames_to_column(as.data.frame(counts(dds, normalized = TRUE)), var = "geneId")
readr::write_tsv(x = normCounts, path = paste0(c(outPrefix,".normCounts.tab"), collapse = ""))

## r-log normalized counts
rld <- rlog(dds, blind = FALSE)
rldCount <- rownames_to_column(as.data.frame(assay(rld)), var = "geneId")
readr::write_tsv(x = rldCount, path = paste0(c(outPrefix,".rlogCounts.tab"), collapse = ""))

###########################################################################

pt_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.key.size = unit(1.2,"cm"),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")
  )

plotPCA(rld, intgroup=c("treatment", "time"), ntop = 4000)

pcaData <- plotPCA(rld, intgroup = c("treatment", "time", "condition"),
                   returnData = TRUE, ntop = 4000)
percentVar <- sprintf("%.2f", 100 * attr(pcaData, "percentVar"))

pcaData$treatment <- forcats::as_factor(pcaData$treatment)
pcaData$time <- forcats::as_factor(pcaData$time)
pcaData$condition <- forcats::as_factor(pcaData$condition)

pltTitle <- "Principal Component Analysis"

fillColumn <- "treatment"
shapeColumn <- "time"

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


pt_pca <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = !!sym(fillColumn), shape = !!sym(shapeColumn)),
             size = 6, stroke = 1) +
  scale_shape_manual(values = c("12hr" = 16, "24hr" = 17)) +
  # guides(fill=guide_legend(override.aes=list(shape=21))) +
  scale_color_manual(values = pointCol) +
  geom_text_repel(mapping = aes(label = name), size = 4,  
                  point.padding = unit(0.5, "lines")) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    title = pltTitle,
    x = paste0("PC1: ",percentVar[1],"% variance"),
    y = paste0("PC2: ",percentVar[2],"% variance")
  ) +
  pt_theme


ggsave(filename = paste(outPrefix, ".PCA.png", sep = ""), device = "png",
       plot = pt_pca, width = 10, height = 8, dpi = 400)


###########################################################################
## PCA based on rld counts

normCountMat <- as.matrix(rldCount[, c(exptInfo$sampleId), drop = FALSE])
rownames(normCountMat) <- rldCount$geneId

## select top 4000 rows with highest variation
rv <- matrixStats::rowVars(normCountMat)
keep <- order(rv, decreasing=TRUE)[seq_len(min(4000, length(rv)))]

## remove low count rows
# keep<- rowSums(normCountMat > 1) >= 2

normCountMatFiltered <- normCountMat[keep, ]

## transform the data such that the genes are columns and each sample is a row
## also append the additional information for each sample using left_join()
df2 <- as.data.frame(t(normCountMatFiltered)) %>% 
  tibble::rownames_to_column(var = "sampleId") %>% 
  dplyr::left_join(y = exptInfo, by = c("sampleId" = "sampleId")) %>% 
  dplyr::select(!!!colnames(exptInfo), dplyr::everything())

row.names(df2) <- df2$sampleId

res.pca <- PCA(df2, graph = FALSE, scale.unit = FALSE,
               quali.sup = 1:ncol(exptInfo), ncp = 10)

eig.val <- get_eigenvalue(res.pca)

## scree plot: variance by PC
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

## Graph of individuals
ind <- get_pca_ind(res.pca)

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
plotData <- as.data.frame(ind$coord) %>%
  tibble::rownames_to_column(var = "sampleId") %>%
  dplyr::left_join(y = exptInfo, by = c("sampleId" = "sampleId"))

pairs(x = plotData[, 2:6],
      pch = 19,  cex = 1,
      col = plotData$condition,
      lower.panel=NULL)


## set the factor levels
plotData$condition <- forcats::as_factor(plotData$condition)
plotData$time <- forcats::as_factor(plotData$time)
plotData$treatment <- forcats::as_factor(plotData$treatment)


pltTitle <- "Principal Component Analysis: rld"

## decide which PCs to use for plotting
pcToPlot <- c(1, 2)
pcCols <- grep(pattern = "Dim.", x = colnames(plotData), value = T)[pcToPlot]
fillColumn <- "treatment"
shapeColumn <- "time"

if(length(unique(plotData[[fillColumn]])) <= 9){
  pointCol <- base::structure(
    .Data = RColorBrewer::brewer.pal(n = length(unique(plotData[[fillColumn]])), name = "Set1"),
    names = levels(plotData[[fillColumn]])
  )
} else{
  pointCol <- base::structure(
    .Data = rainbow(n = length(unique(plotData[[fillColumn]]))),
    names = levels(plotData[[fillColumn]])
  )
}

pointCol <- c("DMSO" = "black", "Tet_5uM" = "#a6cee3", "Tet_15uM" = "#1f78b4",
              "ADTet_5uM" = "#fb9a99", "ADTet_15uM" = "#e31a1c")
shapeCode <- c("12hr" = 16, "24hr" = 17)

pt_rldPca <- ggplot(
  data = plotData,
  mapping = aes(x = !!sym(pcCols[1]), y = !!sym(pcCols[2]), label = sampleId)
) +
  geom_point(mapping = aes(color = !!sym(fillColumn), shape = !!sym(shapeColumn)),
             size = 6, stroke = 2) +
  scale_shape_manual(values = shapeCode) +
  # guides(color=guide_legend(override.aes=list(shape=16))) +
  scale_color_manual(values = pointCol) +
  geom_text_repel(size = 4, point.padding = unit(0.5, "lines")) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) + 
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
  labs(
    title = pltTitle,
    x = paste("PC",pcToPlot[1]," (", sprintf("%.2f", eig.val[pcToPlot[1], "variance.percent"]), "%)", sep = ""),
    y = paste("PC",pcToPlot[2]," (", sprintf("%.2f", eig.val[pcToPlot[2], "variance.percent"]), "%)", sep = "")
  ) +
  pt_theme


ggsave(filename = paste(outPrefix, ".rld_PCA.pdf", sep = ""), device = "pdf",
       plot = pt_rldPca, width = 10, height = 8, dpi = 400)

ggsave(filename = paste(outPrefix, ".rld_PCA.png", sep = ""), device = "png",
       plot = pt_rldPca, width = 10, height = 8, dpi = 400)

## PC scatter plot
pt_pcScatter <- GGally::ggpairs(
  data = plotData, columns = 2:6,
  mapping = aes(color = !!sym(fillColumn), shape = !!sym(shapeColumn)),
  lower = list(continuous = wrap("points", size = 3)),
  upper = "blank",
  diag = list(continuous = wrap("diagAxis", colour = "black")),
  labeller = as_labeller(
    x = structure(
      paste(rownames(eig.val)," (", sprintf("%.2f", eig.val[, "variance.percent"]), "%)", sep = ""),
      names = rownames(eig.val)
    )
  ),
  legend = 6,
  showStrips = FALSE
) +
  scale_color_manual(values = pointCol) +
  labs(title = pltTitle) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.title = element_text(face = "bold", size = 15),
    legend.key.size = unit(1.2,"cm"),
    axis.text = element_blank(),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    plot.title = element_text(size = 16, face = "bold")
  )


ggsave(filename = paste(outPrefix, ".PCA_scatter.png", sep = ""), device = "png",
       plot = pt_pcScatter, width = 12, height = 10, dpi = 400)

###########################################################################
## sample distance matrix
scaledMat <- chipmine::scale_matrix_rows(normCountMatFiltered)
M <- cor(scaledMat)

png(filename = paste(outPrefix, ".zscore_corrplot.png", sep = ""),
    width = 3000, height = 3000, res = 280)
corrplot(corr = M, method = "circle", order = "AOE")
dev.off()

sampleDists <- dist(t(assay(rld)))
# sampleDists <- dist(t(scaledMat))
sampleDistMatrix <- as.matrix(sampleDists)


pt_dist <- ComplexHeatmap::Heatmap(
  matrix = sampleDistMatrix,
  col = colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255),
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 14),
  row_dend_reorder = TRUE, column_dend_reorder = TRUE,
  heatmap_legend_param = list(
    title = "Distance", title_gp = gpar(fontsize = 12),
    title_position = "topcenter")
)


pt_dist <- draw(
  pt_dist,
  column_title = paste("Distance matrix of normalized read counts:", analysisName),
  column_title_gp = gpar(fontface = "bold", fontsize = 16),
  padding = unit(rep(0.5, 4), "cm")
)

png(filename = paste(outPrefix, ".distance_heatmap.png", sep = ""),
    width = 3000, height = 3000, res = 400)
draw(pt_dist)
dev.off()

pdf(file = paste(outPrefix, ".distance_heatmap.pdf", sep = ""),
    width = 8, height = 8)
draw(pt_dist)
dev.off()

###########################################################################
## Poisson dissimilarity matrix: remember it does its own normalization of counts
poisd <- PoiClaClu::PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )

rownames(samplePoisDistMatrix) <- dds$sampleId
colnames(samplePoisDistMatrix) <- dds$sampleId

## use Poisson distance to generate cluster using hclust
poisDistClust <- hclust(poisd$dd)

anDf <- tibble::tibble(sampleId = rownames(samplePoisDistMatrix)) %>% 
  dplyr::left_join(y = exptInfo, by = "sampleId") %>% 
  as.data.frame()

htAn <- ComplexHeatmap::HeatmapAnnotation(
  treatment = anDf$treatment,
  time = anno_simple(
    x = anDf$time, pch = shapeCode[anDf$time],
    col = structure(c("white", "white"), names = names(shapeCode))
  ),
  which = "column",
  # df = dplyr::select(anDf, treatment, time),
  col = list(
    treatment = pointCol
  )
)

## plot heatmap
pt_poisDist <- ComplexHeatmap::Heatmap(
  matrix = samplePoisDistMatrix,
  col = colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(255),
  top_annotation = htAn,
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 14),
  cluster_rows = poisDistClust, cluster_columns = poisDistClust,
  row_dend_reorder = F, column_dend_reorder = F,
  heatmap_legend_param = list(
    title = "Distance", title_gp = gpar(fontsize = 12),
    legend_height = unit(5, "cm"), title_position = "topcenter")
)


pt_poisDist <- draw(
  pt_poisDist,
  column_title = paste("Poisson distance matrix read counts:", analysisName),
  column_title_gp = gpar(fontface = "bold", fontsize = 16),
  merge_legends = TRUE,
  padding = unit(rep(0.5, 4), "cm")
)

png(filename = paste(outPrefix, ".poisson_distance_heatmap.png", sep = ""),
    width = 3500, height = 3000, res = 380)
draw(pt_poisDist)
dev.off()

pdf(file = paste(outPrefix, ".poisson_distance_heatmap.pdf", sep = ""),
    width = 9, height = 8)
draw(pt_poisDist)
dev.off()


#############################################################################
## correlation scatter plot
pt_pairs <- GGally::ggpairs(
  data = dplyr::select(rldCount, -geneId),
  upper = list(continuous = wrap("points", size = 0.1)),
  lower = list(continuous = wrap("cor", size = 10)),
  diag = list(continuous = "densityDiag"),
  title = paste("scatter plot of rlog transformed normalized read counts:", analysisName)
) +
  theme_bw() +
  theme(
    strip.text.y = element_text(size = 18, angle = 0, hjust = 0, face = "bold"),
    strip.text.x = element_text(size = 18, angle = 90, hjust = 0, face = "bold")
  )

png(filename = paste(outPrefix, ".scatter_matrix.png", sep = ""),
    width = 10000, height = 10000, res = 300)

print(pt_pairs)
dev.off()

# ###########################################################################
# ## calculate mean FPKM
# 
# exptInfo <- read.table(file = sampleInfoFile, header = T, sep = "\t", row.names = "sampleId")
# 
# ## IF NEEDED: modify the row names and column names
# # colnames(rawData) <- sub("_WT", "_", colnames(rawData))
# # rownames(exptInfo) <- sub("_WT", "_", rownames(exptInfo))
# 
# ## select only those sample rows which are part of current comparison
# designInfo <- exptInfo
# # designInfo <- droplevels(subset(exptInfo, condition %in% compare))
# 
# if(! all( rownames(designInfo) %in% colnames(geneExpression) ) ){
#   stop("Column names in FPKM matrix does not match with row names in experiment data")
# }
# 
# 
# ## get the sample names for each condition under comparison
# grp1 <- sapply(rownames(designInfo[designInfo$condition %in% compare[1], ]), FUN = as.name, USE.NAMES = F, simplify = T)
# name1 <- paste(compare[1], "_meanFPKM", sep = "")
# grp1Len <- length(grp1)
# 
# grp2 <- sapply(rownames(designInfo[designInfo$condition %in% compare[2], ]), FUN = as.name, USE.NAMES = F, simplify = T)
# name2 <- paste(compare[2], "_meanFPKM", sep = "")
# grp2Len <- length(grp2)
# 
# ## can add multiple groups
# 
# ## calculate the average FPKM
# fpkmData <- geneExpression %>% 
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



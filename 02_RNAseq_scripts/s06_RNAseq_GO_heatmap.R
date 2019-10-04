library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)


## volcano plot of RNAseq DEGs with specific GO term genes marked
## heatmap of a geneset belonging to GO Terms


rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")
source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/topGO_functions.R")

###########################################################################

analysisName <- "PG_HE_8mpf_vs_PG_WT_8mpf"

degResults <- c(analysisName)

outDir <- here::here("analysis", "02_DESeq2_diff", analysisName, "enrichment_plots")

outPrefix <- paste(outDir, "/", analysisName, ".GO", sep = "")

file_sampleInfo <- here::here("data", "sample_info.txt")
file_RNAseq_info <- here::here("data", "RNAseq_info.txt", sep = "")
diffDataPath <- here::here("analysis", "02_DESeq2_diff")

file_termSet <- paste(outDir, "/", "go_heatmap_terms.txt", sep = "")

lfc_col = "shrinkLog2FC"
fdr_col = "padj"
###########################################################################
## prepare data
sampleInfo <- suppressMessages(readr::read_tsv(file = file_sampleInfo))

termSet <- suppressMessages(readr::read_tsv(file = file_termSet))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == degResults)

sampleIds <- unlist(strsplit(x = rnaseqInfo$samples, split = ","))

degs <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[1])) 

goRes <- suppressMessages(readr::read_tsv(file = rnaseqInfo$topGO[1])) %>% 
  dplyr::filter(GO.ID %in% termSet$id) %>% 
  tidyr::separate_rows(genes, sep = ",") %>% 
  dplyr::select(GO.ID, Term, genes)

goRes$Term <- stringr::str_wrap(str = goRes$Term, width = 30)

degs <- dplyr::left_join(x = degs, y = goRes, by = c("geneId" = "genes")) %>% 
  tidyr::replace_na(replace = list(Term = "none"))

normCounts <- suppressMessages(readr::read_tsv(file = rnaseqInfo$normCount[1])) %>% 
  dplyr::select(geneId, sampleIds)

rldCounts <- suppressMessages(readr::read_tsv(file = rnaseqInfo$rld[1])) %>% 
  dplyr::select(geneId, sampleIds)

###########################################################################
## volcano plot
vp <- volcano_plot(df = degs, fdr_col = fdr_col, lfc_col = lfc_col,
                   ylimit = 40, xlimit = c(-5, 5))


termColor <- structure(
  .Data = RColorBrewer::brewer.pal(n = length(unique(goRes$Term)), name = "Set1"), 
  names = unique(goRes$Term)
)


vpt <- ggplot(mapping = aes(x = !!sym(lfc_col), y = log10FDR)) +
  geom_point(data = vp$data, color = "grey", alpha = 0.5) +
  geom_point(data = dplyr::filter(vp$data, Term != "none"),
             mapping = aes(color = Term),
             size = 2, alpha = 0.7) +
  # geom_text_repel(data = dplyr::filter(vp$data, Term != "none"),
  #                  mapping = aes(label = GENE_NAME)) +
  scale_color_manual(
    name = "GO Term",
    values = termColor
  ) +
  scale_x_continuous(expand = expand_scale(add = 0.1)) +
  scale_y_continuous(expand = expand_scale(add = 1)) +
  labs(
    title = paste("Volcano plot:", analysisName),
    x = "log2(fold-change)",
    y = "-log10(q-value)"
  ) +
  theme_bw() +
  theme(legend.text = element_text(size = 12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_text(face = "bold", size = 14),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 14))

pdf(file = paste(outPrefix, ".volcano.pdf", sep = ""), width = 10, height = 10)
vpt
dev.off()


###########################################################################
## heatmap

geneSubset <- dplyr::filter(degs, Term != "none") %>% 
  dplyr::left_join(y = rldCounts, by = "geneId") %>% 
  dplyr::mutate(id = row_number())

rownameCol <- "GENE_NAME"
showRowNames <- TRUE
rowNameFontSize <- 14
colNameFontSize <- 14


## fold change heatmap
foldChangeDf <- dplyr::select(geneSubset, id, !!lfc_col) %>%
  tibble::column_to_rownames(var = "id")

foldChangeMat <- data.matrix(foldChangeDf)

fcHeatmap <- Heatmap(
  matrix = foldChangeMat,
  col = colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"), space = "LAB"),
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_labels = geneSubset$GENE_NAME,
  row_names_gp = gpar(fontsize = rowNameFontSize),
  column_names_gp = gpar(fontsize = colNameFontSize), 
  width = unit(2, "cm"),
  heatmap_legend_param = list(title = "\nlog2(fold_change)")
)


geneCounts <- dplyr::select(geneSubset, id, !!!sampleIds) %>%
  column_to_rownames(var = "id")

countMat <- data.matrix(geneCounts)

countZscoreMat <- chipmine::scale_matrix_rows(x = countMat)


## plot main rld score heatmap
rldZscoreHeatmap <- Heatmap(
  matrix = countZscoreMat,
  name = "countHt",
  col = colorRamp2(breaks = c(min(countZscoreMat), 0, max(countZscoreMat)),
                   colors = c("green", "black", "red"), space = "LAB"), 
  row_names_gp = gpar(fontsize = rowNameFontSize),
  column_names_gp = gpar(fontsize = colNameFontSize), 
  cluster_columns = FALSE,
  row_title_rot = 0,
  width = unit(10, "cm"), row_names_max_width = unit(15, "cm"), 
  heatmap_legend_param = list(title = "z-score(rlog counts)", color_bar = "continuous")
) 


htAn <- rowAnnotation(
  goTerm = geneSubset$Term,
  col = list(goTerm = termColor),
  show_legend = FALSE
)

## log2(fold_change) + rld heatmap with annotation
htList <- htAn + rldZscoreHeatmap + fcHeatmap 


# png(filename = paste(outPrefix, "_fc_rld_heatmap.png", sep = ""), width=6000, height=6000, res = 550)
pdf(file = paste(outPrefix, ".heatmap.pdf", sep = ""), width = 12, height = 14)
draw(object = htList,
     main_heatmap = "countHt",
     column_title = paste("Heatmap of DEGs belonging to GO terms:", analysisName),
     split = geneSubset$Term,
     show_row_dend = FALSE,
     # auto_adjust = FALSE,
     row_sub_title_side = "left",
     column_title_gp = gpar(fontsize = 14)
)

dev.off()





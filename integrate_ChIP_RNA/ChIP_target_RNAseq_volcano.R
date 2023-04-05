library(chipmine)
library(org.AFumigatus.Af293.eg.db)
library(ggpubr)
library(ggrepel)
library(cowplot)


rm(list = ls())
# source("E:/Chris_UM/GitHub/omics_util/RNAseq_scripts/DESeq2_functions.R")
source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")


analysisName <- "5A9_C_vs_CEA17_C"
diffPair <- "5A9_C_vs_CEA17_C"
chipSamples <- c("CREEHA_CONTROL4", "CREEHA_10MMAA4")

file_ChIPtargets <- here::here("analysis", "ChIPseq_analysis", "peak_targets", "peak_targets.curated.filtered.tab")
file_degs <- here::here("analysis", "RNAseq_data", diffPair, "5A9_C_vs_CEA17_C.DEG_all.txt")

outDir <- here::here("analysis", "integration_analysis", "ChIPseq_RNAseq_volcano", analysisName)

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, analysisName, sep = "/")

orgDb <- org.AFumigatus.Af293.eg.db
file_goMap <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_orgDb/geneid2go.AFumigatus_Af293.topGO.map"

##################################################################################

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered"),
  FUN = function(x){ structure(paste(x, ".", chipSamples, sep = ""), names = chipSamples) },
  simplify = F, USE.NAMES = T
)


FDR_cut <- 0.05
lfc_cut <- 0.585
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

##################################################################################
## prepare dataset
targetSet <- suppressMessages(readr::read_tsv(file = file_ChIPtargets, col_names = T))

diffData <- suppressMessages(readr::read_tsv(file = file_degs, col_names = T))


fdr_col = "padj"
lfc_col = "shrinkLog2FC"
ylimit = 50
xlimit = c(-5, 5)


plotData <- dplyr::left_join(x = diffData, y = targetSet, by = c("geneId" = "geneId")) %>% 
  tidyr::replace_na(purrr::set_names(list(FALSE, FALSE), unname(tfCols$hasPeak[chipSamples]))) %>% 
  dplyr::mutate(
    category = dplyr::case_when(
      !! as.name(fdr_col) < !! FDR_cut & !! as.name(lfc_col) >= !! up_cut ~ "Significant Up",
      !! as.name(fdr_col) < !! FDR_cut & !! as.name(lfc_col) <= !! down_cut ~ "Significant Down",
      !! as.name(fdr_col) < !! FDR_cut ~ "Significant",
      TRUE ~ "Non-significant"
    )
  )

## store data
readr::write_tsv(x = plotData, path = paste(outPrefix, "_ChIPseq_targets.tab"))


##################################################################################
## generate summary stats
summary_tf1 <- dplyr::filter(plotData, !! as.name(unname(tfCols$hasPeak[chipSamples[1]])) == TRUE) %>% 
  dplyr::group_by(category) %>% 
  summarise(!! unname(tfCols$hasPeak[chipSamples[1]]) := n()) %>% 
  dplyr::ungroup()

summary_tf2 <- dplyr::filter(plotData, !! as.name(unname(tfCols$hasPeak[chipSamples[2]])) == TRUE) %>% 
  dplyr::group_by(category) %>% 
  summarise(!! unname(tfCols$hasPeak[chipSamples[2]]) := n()) %>% 
  dplyr::ungroup()

summary_tf1_specific <- dplyr::filter(
  plotData,
  !! as.name(unname(tfCols$hasPeak[chipSamples[1]])) == TRUE,
  !! as.name(unname(tfCols$hasPeak[chipSamples[2]])) == FALSE) %>% 
  dplyr::group_by(category) %>% 
  summarise(!! chipSamples[1] := n()) %>% 
  dplyr::ungroup()

summary_tf2_specific <- dplyr::filter(
  plotData,
  !! as.name(unname(tfCols$hasPeak[chipSamples[2]])) == TRUE,
  !! as.name(unname(tfCols$hasPeak[chipSamples[1]])) == FALSE) %>% 
  dplyr::group_by(category) %>% 
  summarise(!! chipSamples[2] := n()) %>% 
  dplyr::ungroup()

summary_common <- dplyr::filter_at(plotData, .vars = vars(starts_with("hasPeak.")),
                                   .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::group_by(category) %>% 
  summarise(common = n()) %>% 
  dplyr::ungroup()

summary_combined <- dplyr::group_by(plotData, category) %>% 
  summarise(all_genes = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(y = summary_tf1_specific, by = "category") %>% 
  dplyr::left_join(y = summary_common, by = "category") %>% 
  dplyr::left_join(y = summary_tf2_specific, by = "category")

summary_table <- dplyr::bind_rows(
  summary_combined,
  dplyr::summarise_if(.tbl = summary_combined, .predicate = is.numeric, .funs = sum, na.rm = TRUE) %>% 
    mutate(category = "total")) %>% 
  ggpubr::ggtexttable(rows = NULL, theme = ttheme(base_style = "classic", base_size = 10))


##################################################################################
## volcano plot
plotTitle <- paste(diffPair, "DEGs marked with ChIPseq target genes", sep = " ")

markGenes <- plotData$geneId[ which(plotData[, tfCols$hasPeak[chipSamples[1]], drop = TRUE ] ) ]

## squish the value to limits
plotData$log10FDR <- -log10(plotData[[fdr_col]])
plotData$log10FDR <- scales::squish(x = plotData$log10FDR, range = c(0, ylimit))
plotData[[lfc_col]] <- scales::squish(x = plotData[[lfc_col]], range = xlimit)

## show any specific category only
significantData <- dplyr::filter(plotData, category %in% c("Significant Up", "Significant Down"))

#draw Volcano plot
vpt <- ggplot(mapping = aes(x = !! as.name(lfc_col), y = log10FDR)) +
  geom_hline(yintercept = -log10(FDR_cut), color = "black", linetype = "dashed") +
  geom_vline(xintercept = -lfc_cut, color = "black", linetype = "dashed") +
  geom_vline(xintercept = lfc_cut, color = "black", linetype = "dashed") +
  geom_point(data = plotData, color = "grey", alpha=0.5, size=1.5) +
  geom_point(
    data = dplyr::filter(significantData, !! as.name(tfCols$hasPeak[chipSamples[1]]) == TRUE),
    mapping = aes(color = !!chipSamples[1]), alpha=0.7, size=1.5
  ) +
  geom_point(
    data = dplyr::filter(significantData, !! as.name(tfCols$hasPeak[chipSamples[2]]) == TRUE),
    mapping = aes(color = !!chipSamples[2]), alpha=0.7, size=1.5
  ) +
  geom_point(
    data = dplyr::filter_at(significantData, vars(unname(tfCols$hasPeak[chipSamples])), all_vars(. == TRUE)),
    mapping = aes(color = "common"), size=1.5
  ) +
  scale_color_manual(
    name = "ChIPseq targets",
    values = structure(c("green", "blue", "red"), names = c(chipSamples, "common")),
    breaks = c(chipSamples, "common")
    
  ) +
  scale_x_continuous(name = "log2(fold_change)", limits = xlimit, expand = expand_scale(mult = 0.02)) +
  scale_y_continuous(name = "-log10(q-value)", limits = c(0, ylimit), expand = expand_scale(mult = 0.02)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(legend.background = element_rect(colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.margin = unit(rep(0.5, 4),"cm"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14)) +
  ggtitle(plotTitle)


## mark gene of interest
markGenes <- c()
showNames <- TRUE

## highlight genes of interest
if(length(markGenes) > 0){
  tmpDf <- dplyr::filter(plotData, geneId %in% markGenes)
  
  ## draw the points
  vpt <- vpt + geom_point(data = tmpDf, color = "black", shape = 1)
  
  ## show the gene lables
  if(isTRUE(showNames)){
    vpt <- vpt +
      geom_text_repel(data = tmpDf, mapping = aes(label = geneId),
                      segment.color = '#cccccc',
                      segment.size = 1,
                      size = 4) 
  }
}


mergedPt <- ggarrange(
  vpt, summary_table,
  nrow = 2, heights = c(1, 0.2)
)


png(filename = paste(outPrefix, "_ChIPseq_targets.significant.volcano.png", sep = ""), width = 3000, height = 3000, res = 320)
plot(mergedPt)
dev.off()

pdf(file = paste(outPrefix, "_ChIPseq_targets.significant.volcano.pdf", sep = ""), width = 10, height = 10)
plot(mergedPt)
dev.off()

##################################################################################
## GO analysis: DEG down + peak targets
nodeSize <- 5

plotTitle <- paste("GO enrichment: DEG-Down", diffPair, " and peak targets")
downDEGTargets <- dplyr::filter_at(.tbl = significantData,
                                   .vars = vars(starts_with("hasPeak.")),
                                   .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::filter(category == "Significant Down")


goDataDown <- topGO_enrichment(genes = downDEGTargets$geneId, goMapFile = file_goMap, goNodeSize = nodeSize)

readr::write_tsv(x = goDataDown, path = paste(outPrefix, ".ChIPseq_targets.down.topGO.tab", sep = ""))

topGoScatterDown <- NULL
if(nrow(goDataDown) > 0){
  topGoScatterDown <- topGO_scatterPlot(df = goDataDown, title = plotTitle)
}

plotOut <- paste(outPrefix, ".ChIPseq_targets.down.topGO.png", sep = "")

# topGO <- go_and_scatterPlot(genes = downDEGTargets$geneId,
#                             goToGeneFile = file_goMap,
#                             goTitle = plotTitle,
#                             plotOut = plotOut, goNodeSize = nodeSize)

##################################################################################
## GO analysis: DEG up + peak targets
plotTitle <- paste("GO enrichment: DEG-Up", diffPair, " and peak targets")
upDEGTargets <- dplyr::filter_at(.tbl = significantData,
                                   .vars = vars(starts_with("hasPeak.")),
                                   .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::filter(category == "Significant Up")


goDataUp <- topGO_enrichment(genes = upDEGTargets$geneId, goMapFile = file_goMap, goNodeSize = nodeSize)

readr::write_tsv(x = goDataUp, path = paste(outPrefix, ".ChIPseq_targets.up.topGO.tab", sep = ""))

topGoScatterUp <- NULL
if(nrow(goDataUp) > 0){
  topGoScatterUp <- topGO_scatterPlot(df = goDataUp, title = plotTitle)
}

plotOut <- paste(outPrefix, ".ChIPseq_targets.down.up.png", sep = "")

# topGO <- go_and_scatterPlot(genes = upDEGTargets$geneId,
#                             goToGeneFile = file_goMap,
#                             goTitle = plotTitle,
#                             plotOut = plotOut, goNodeSize = nodeSize)


##################################################################################
aligned_plots <- align_plots(
  plotlist = list(down = topGoScatterDown, up = topGoScatterUp),
  align = "v")

pdf(file = paste(outPrefix, ".ChIPseq_targets.topGO.pdf", sep = ""), width = 12, height = 12, onefile = T, pointsize = 18)
ggdraw(aligned_plots$down)
ggdraw(aligned_plots$up)
dev.off()






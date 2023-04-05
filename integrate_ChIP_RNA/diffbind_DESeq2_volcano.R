library(chipmine)
library(org.AFumigatus.Af293.eg.db)
library(ggpubr)
library(ggrepel)
library(cowplot)


rm(list = ls())
# source("E:/Chris_UM/GitHub/omics_util/RNAseq_scripts/DESeq2_functions.R")
# source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")

# "CEA17_AA_vs_CEA17_C", "5A9_AA_vs_5A9_C", "5A9_C_vs_CEA17_C", "5A9_AA_vs_CEA17_AA"
analysisName <- "5A9_AA_vs_5A9_C"
diffPair <- "5A9_AA_vs_5A9_C"

outDir <- here::here("analysis", "integration_analysis", "diffbind_DESeq2_volcano", analysisName)
outPrefix <- paste(outDir, "/", analysisName, ".diffbind_goodPeaks", sep = "")

##################################################################################
chipSamples <- c("CREEHA_CONTROL4", "CREEHA_10MMAA4")
diffbindCompare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

file_diffbindTargets <- here::here("analysis", "ChIPseq_analysis",
                                   "peak_targets", "diffbind_allPeak_targets.tab")
file_deseq2 <- paste(diffPair, ".DEG_all.txt", sep = "")
file_degs <- here::here("analysis", "RNAseq_data", diffPair, file_deseq2)


if(!dir.exists(outDir)){
  dir.create(path = outDir, recursive = TRUE)
}


orgDb <- org.AFumigatus.Af293.eg.db
file_goMap <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_orgDb/geneid2go.AFumigatus_Af293.topGO.map"

##################################################################################

grp1 <- diffbindCompare[1]
grp2 <- diffbindCompare[2]
grp1Enrich = paste(grp1, ":enriched", sep = "")
grp2Enrich = paste(grp2, ":enriched", sep = "")
grp1Specific = paste(grp1, ":specific", sep = "")
grp2Specific = paste(grp2, ":specific", sep = "")

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered"),
  FUN = function(x){ structure(paste(x, ".", chipSamples, sep = ""), names = chipSamples) },
  simplify = F, USE.NAMES = T
)

fdr_col = "padj"
lfc_col = "shrinkLog2FC"
ylimit = 50
xlimit = c(-5, 5)

FDR_cut <- 0.05
lfc_cut <- 0.585
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

##################################################################################
## prepare dataset
diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffbindTargets, col_names = T)) %>% 
  dplyr::select(-starts_with("summitSeq."))

diffbindRes <- dplyr::filter(diffbindRes, pvalFilteredN > 0)

degData <- suppressMessages(readr::read_tsv(file = file_degs, col_names = T))


mergedData <- dplyr::left_join(x = degData, y = diffbindRes, by = c("geneId" = "geneId")) %>% 
  tidyr::replace_na(purrr::set_names(list(FALSE, FALSE), unname(tfCols$hasPeak[chipSamples]))) %>% 
  dplyr::mutate(
    categoryRNAseq = dplyr::case_when(
      !! as.name(fdr_col) < !! FDR_cut & !! as.name(lfc_col) >= !! up_cut ~ "Significant Up",
      !! as.name(fdr_col) < !! FDR_cut & !! as.name(lfc_col) <= !! down_cut ~ "Significant Down",
      !! as.name(fdr_col) < !! FDR_cut ~ "Significant",
      TRUE ~ "Non-significant"
    )
  )

## store data
readr::write_tsv(x = mergedData, path = paste(outPrefix, ".data.tab"))

##################################################################################
## select best peak for a gene which has multiple peaks
plotData <- mergedData %>% 
  dplyr::mutate(
    minPeakDist = pmin(abs(!!as.name(tfCols$peakDist[chipSamples[1]])),
                       abs(!!as.name(tfCols$peakDist[chipSamples[2]])),
                       na.rm = TRUE
                       )
  ) %>% 
  dplyr::group_by(geneId) %>% 
  dplyr::arrange(abs(minPeakDist), desc(bestPval)) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup()

## generate summary stats
summary_tf1 <- dplyr::filter(plotData, !! as.name(unname(tfCols$hasPeak[chipSamples[1]])) == TRUE) %>% 
  dplyr::group_by(categoryRNAseq) %>% 
  summarise(!! unname(tfCols$hasPeak[chipSamples[1]]) := n()) %>% 
  dplyr::ungroup()

summary_tf2 <- dplyr::filter(plotData, !! as.name(unname(tfCols$hasPeak[chipSamples[2]])) == TRUE) %>% 
  dplyr::group_by(categoryRNAseq) %>% 
  summarise(!! unname(tfCols$hasPeak[chipSamples[2]]) := n()) %>% 
  dplyr::ungroup()

summary_tf1_enrich <- dplyr::filter(plotData, categoryDiffbind %in% c(grp1Enrich, grp1Specific)) %>% 
  dplyr::group_by(categoryRNAseq) %>% 
  summarise(!! chipSamples[1] := n()) %>% 
  dplyr::ungroup()

summary_tf2_enrich <- dplyr::filter(plotData, categoryDiffbind %in% c(grp2Enrich, grp2Specific)) %>% 
  dplyr::group_by(categoryRNAseq) %>% 
  summarise(!! chipSamples[2] := n()) %>% 
  dplyr::ungroup()

summary_common <- dplyr::filter(plotData, categoryDiffbind == "common") %>% 
  dplyr::group_by(categoryRNAseq) %>% 
  summarise(common = n()) %>% 
  dplyr::ungroup()

summary_combined <- dplyr::group_by(plotData, categoryRNAseq) %>% 
  summarise(all_genes = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(y = summary_tf1_enrich, by = "categoryRNAseq") %>% 
  dplyr::left_join(y = summary_common, by = "categoryRNAseq") %>% 
  dplyr::left_join(y = summary_tf2_enrich, by = "categoryRNAseq")

summary_table <- dplyr::bind_rows(
  summary_combined,
  dplyr::summarise_if(.tbl = summary_combined, .predicate = is.numeric, .funs = sum, na.rm = TRUE) %>% 
    mutate(categoryRNAseq = "total")) %>% 
  ggpubr::ggtexttable(rows = NULL, theme = ttheme(base_style = "classic", base_size = 10))



##################################################################################


plotTitle <- paste(diffPair, "DEGs marked with ChIPseq target DiffBind change", sep = " ")

markGenes <- plotData$geneId[ which(plotData[, tfCols$hasPeak[chipSamples[1]], drop = TRUE ] ) ]

## squish the value to limits
plotData$log10FDR <- -log10(plotData[[fdr_col]])
plotData$log10FDR <- scales::squish(x = plotData$log10FDR, range = c(0, ylimit))
plotData[[lfc_col]] <- scales::squish(x = plotData[[lfc_col]], range = xlimit)

## show any specific categoryRNAseq only
significantData <- dplyr::filter(plotData, categoryRNAseq %in% c("Significant Up", "Significant Down"))
# significantData <- plotData

# targetColors <- c("specific:CREEHA_CONTROL:down" = "#a50026",
#                   "specific:CREEHA_CONTROL:noDiff" = "#f46d43",
#                   "common:down" = "#fee08b",
#                   "common:noDiff" = "yellow",
#                   "common:up" = "#e0f3f8",
#                   "specific:CREEHA_10MMAA:noDiff" = "#74add1",
#                   "specific:CREEHA_10MMAA:up" = "#006837")
# 
# diffbindColors <- c("down" = "red", "noDiff" = "green", "up" = "blue")
# targetTypeShape <- c("specific:CREEHA_CONTROL" = 25,
#                      "common" = 22,
#                      "specific:CREEHA_10MMAA" = 24)

diffbindCategoryColors <- structure(
  c("#e60000", "#ff4d4d", "green", "#6666ff", "#0000e6"),
  names = c(grp1Specific, grp1Enrich, "common", grp2Enrich, grp2Specific))



#draw Volcano plot
vpt <- ggplot(mapping = aes(x = !! as.name(lfc_col), y = log10FDR)) +
  geom_hline(yintercept = -log10(FDR_cut), color = "black", linetype = "dashed") +
  geom_vline(xintercept = -lfc_cut, color = "black", linetype = "dashed") +
  geom_vline(xintercept = lfc_cut, color = "black", linetype = "dashed") +
  geom_point(data = plotData, color = "grey", alpha=0.5, size=2) +
  geom_point(
    data = dplyr::filter(significantData, !is.na(categoryDiffbind)),
    mapping = aes(color = categoryDiffbind), alpha=0.7, size=2
  ) +
  scale_color_manual(
    name = "DiffBind change", values = diffbindCategoryColors,
    breaks = names(diffbindCategoryColors)
  ) +
  # geom_point(
  #   data = dplyr::filter(significantData, !is.na(categoryDiffbind)),
  #   mapping = aes(color = diffBind, fill = diffBind, shape = peakOccupancy),
  #   alpha=0.7, size=2
  # ) +
  # scale_color_manual(
  #   name = "DiffBind change", values = diffbindColors, breaks = names(diffbindColors)
  # ) +
  # scale_fill_manual(
  #   values = diffbindColors, breaks = names(diffbindColors), guide = FALSE
  # ) +
# scale_shape_manual(
#   name = "peak occupancy", values = targetTypeShape, breaks = names(targetTypeShape)
# ) +
scale_x_continuous(name = "log2(fold_change)", limits = xlimit, expand = expand_scale(mult = 0.02)) +
  scale_y_continuous(name = "-log10(q-value)", limits = c(0, ylimit), expand = expand_scale(mult = 0.02)) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5, fill = "black"))) +
  theme_bw() +
  theme(legend.background = element_rect(colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.margin = unit(rep(0.5, 4),"cm"),
        axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 30)) +
  ggtitle(plotTitle)

# plot(vpt)

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


png(filename = paste(outPrefix, ".significant.volcano.png", sep = ""), width = 3500, height = 3500, res = 300)
plot(mergedPt)
dev.off()

pdf(file = paste(outPrefix, ".significant.volcano.pdf", sep = ""), width = 12, height = 12)
plot(mergedPt)
dev.off()




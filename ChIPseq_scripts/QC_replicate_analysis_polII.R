suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))


rm(list = ls())


##################################################################################
analysisName <- "PolII_replicate"
outDir <- here::here("analysis", "02_QC_polII")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_replicates <- here::here("analysis", "02_QC_polII", "polII_replicates.txt")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")

file_genes <- here::here("data", "reference_data", "AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

polII_dataPath <- here::here("data", "polII_data")

file_polII <- paste(polII_dataPath, "/", "samples_polII.all.list", sep = "")

##################################################################################

geneSet <- suppressMessages(readr::read_tsv(
  file = file_genes,
  col_names = c("chr", "start", "end", "geneId", "score", "strand"))) %>%
  dplyr::select(geneId)

polIIsamples <- readr::read_tsv(file = file_polII,  comment = "#")

polIIInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = polIIsamples$sampleId,
  dataPath = polII_dataPath)


replicatePairs <- suppressMessages(readr::read_tsv(file = file_replicates))

plotListAll <- list()
plotList_pval_distibution <- list()

i <- 1

corrDf <- NULL

pdf(file = paste(outPrefix, ".correlation.pdf", sep = ""), width = 15, height = 10,
    onefile = TRUE, pointsize = 8)


for (i in 1:nrow(replicatePairs)) {
  
  cat(i, ":", replicatePairs$rep1[i], "\n")
  
  repInfo <- dplyr::filter(polIIInfo, sampleId %in% c(replicatePairs$rep1[i], replicatePairs$rep2[i]))
  
  if(nrow(repInfo) == 0){
    next
  }
  
  if(nrow(repInfo) != 2){
    stop("sampleInfo should have 2 rows for 2 replicates to be compared")
  }
  
  rep1Col <- repInfo$sampleId[1]
  rep2Col <- repInfo$sampleId[2]
  exprDf <- get_polII_expressions(genesDf = geneSet, exptInfo = repInfo) 
  
  pairCor <- compare_replicates(data = exprDf, rep1Col = rep1Col, rep2Col = rep2Col,
                                trans = "log2", pseudoCount = 0.2)
  
  plot(pairCor$figure)
  
  pairCorrDf <- dplyr::mutate(.data = pairCor$data, rep1 = rep1Col, rep2 = rep2Col)
  corrDf <- dplyr::bind_rows(corrDf, pairCorrDf)
  
  Sys.sleep(5)
}

dev.off()

##################################################################################
corrDf <- dplyr::mutate(corrDf, fractionPer = forcats::as_factor(fractionPer)) %>% 
  dplyr::left_join(
    y = dplyr::select(polIIInfo, sampleId, SM_TF, timePoint), by = c("rep1" = "sampleId")
  )

cummulativeCorr <- tidyr::pivot_wider(
  data = corrDf,
  id_cols = c(rep1, rep2, SM_TF, timePoint),
  names_from = fractionPer,
  values_from = c(pearson, spearman)
) %>% 
  dplyr::arrange(SM_TF, desc(`pearson_100%`)) %>% 
  dplyr::select(SM_TF, timePoint, rep1, rep2, everything())

readr::write_tsv(
  x = cummulativeCorr,
  path = paste(outPrefix, ".cummulative_correlation.tab", sep = "")
)



gg_pearsonScatter <- dplyr::filter(corrDf, fractionPer %in% c("100%", "50%", "10%")) %>% 
  ggplot(mapping = aes(x = fractionPer, y = pearson, fill = pearson)) +
  geom_boxplot(width=0.4, fill = NA, outlier.colour = NA, color = alpha("black", 0.7)) +
  geom_quasirandom(size = 3, shape = 21) +
  scale_fill_distiller(limits = c(-1, 1), type = "seq", palette = "RdBu") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Pearson correlation for top X% genes in polII ChIPseq replicates",
    subtitle = "replicate average signal is used to rank the genes and then top X% genes are selected",
    x = "top X% genes", y = "Pearson correlation"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  )


ggsave(filename = paste(outPrefix, ".pearson_scatter.png", sep = ""),
       plot = gg_pearsonScatter, width = 8, height = 8)
ggsave(filename = paste(outPrefix, ".pearson_scatter.pdf", sep = ""),
       plot = gg_pearsonScatter, width = 8, height = 8)







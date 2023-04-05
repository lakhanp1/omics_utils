suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))


rm(list = ls())

##################################################################################
analysisName <- "TF_replicate"
outDir <- here::here("analysis", "02_QC_TF")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_replicates <- here::here("analysis", "02_QC_TF", "tf_replicates.txt")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")

file_genes <- here::here("data", "reference_data", "AN_genes_for_polII.bed")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")

matrixType <- "2kb_summit"
up <- 2000
down <- 2000
body <- 0
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)

##################################################################################

tfSampleList <- readr::read_tsv(file = file_tf_macs2, comment = "#")

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$sampleId,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType)

# tfInfoList <- purrr::transpose(tfInfo)  %>% 
#   purrr::set_names(nm = purrr::map(., "sampleId"))


replicatePairs <- suppressMessages(readr::read_tsv(file = file_replicates))

plotListAll <- list()
plotList_pval_distibution <- list()
i <- 1

pdf(file = paste(outPrefix, ".correlation.pdf", sep = ""), width = 18, height = 12,
    onefile = TRUE, pointsize = 8)

for (i in 1:nrow(replicatePairs)) {
  # for (i in 1:10) {
  
  cat(i, ":", replicatePairs$rep1[i], "\n")
  
  repInfo <- dplyr::filter(tfInfo, sampleId %in% c(replicatePairs$rep1[i], replicatePairs$rep2[i]))
  if(nrow(repInfo) == 0){
    next
  }
  
  plots_pval <- compare_ChIPseq_replicates(
    sampleInfo = repInfo, compare = "pvalue", yintercept = 20,
    peakFormat = repInfo$peakType[1]
  )
  
  plots_enrichment <- compare_ChIPseq_replicates(
    sampleInfo = repInfo, compare = "enrichment", yintercept = 3,
    peakFormat = repInfo$peakType[1]
  )
  
  ## all plots combined in a row for a replicate
  repPlot <- ggarrange(
    plots_pval$table,
    ggarrange(
      plots_pval$distribution, plots_pval$valueScatter, plots_pval$rankScatter,
      nrow = 1),
    ggarrange(
      plots_enrichment$distribution, plots_enrichment$valueScatter,
      plots_enrichment$rankScatter, nrow = 1),
    nrow = 3, heights = c(0.1, 0.45, 0.45)) +
    theme(plot.background = element_rect(color = "black"),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
  
  plot(repPlot)
  # plotListAll[[i]] <- repPlot
  
  ## summary table and distribution plot combined
  pvalDistrPlot <- ggarrange(
    plots_pval$table, plots_pval$distribution,
    ncol = 1, heights = c(1, 8)
  ) +
    theme(plot.background = element_rect(color = "black"),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
  
  plotList_pval_distibution[[i]] <- pvalDistrPlot
}

dev.off()


fig_pval_distribution <- ggarrange(plotlist = plotList_pval_distibution,
                                   nrow = 8, ncol = 11, hjust = 0.5)

png(filename = paste(outPrefix, ".pval_distribution.png", sep = ""),
    width = 12000, height = 9000, res = 160)
fig_pval_distribution
dev.off()







suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(here))


## get summit sequence for AN0153 ChIPseq peaks

rm(list = ls())

##################################################################################

analysisName <- "AN0153"
outDir <- here("03_ChIPseq_scripts", "testData")
outPrefix <- paste(outDir, "/", analysisName, sep = "")


orgDb <- org.Anidulans.FGSCA4.eg.db
genome <- BSgenome.Anidulans.FGSCA4.AspGD

peakFile <- here("03_ChIPseq_scripts", "testData", "AN0153_sCopy_OE_16h_HA_ChIPMix46_1.withoutCtrl_peaks.narrowPeak")

##################################################################################

summit200 <- get_peak_summit_seq(
  file = peakFile,
  peakFormat = "narrowPeak",
  genome = genome, length = 200,
  column_suffix = FALSE
)

summit200 <- dplyr::select(summit200, peakId, summitSeq) %>% 
  dplyr::rename(summitSeq200 = summitSeq)

summit500 <- get_peak_summit_seq(
  file = peakFile,
  peakFormat = "narrowPeak",
  genome = genome, length = 500,
  column_suffix = FALSE
)

summit500 <- dplyr::select(summit500, peakId, summitSeq) %>% 
  dplyr::rename(summitSeq500 = summitSeq)


peakDf <- import_peaks_as_df(file = peakFile, column_suffix = FALSE)

peakDf <- dplyr::left_join(x = peakDf, y = summit200, by = "peakId") %>% 
  dplyr::left_join(y = summit500, by = "peakId")

readr::write_tsv(x = peakDf, path = paste(outPrefix, "-summit_seq.tab", sep = ""))
















library(chipmine)
library(org.Anidulans.eg.db)
library(TxDb.Anidulans.AspGD.GFF)
library(foreach)
library(doParallel)
library(here)
library(GGally)


rm(list = ls())
# cl <- makeCluster(4) #not to overload your computer
# registerDoParallel(cl)

##################################################################################

file_exptInfo <- here::here("data", "referenceData/sample_info.txt")

file_genes <- here::here("data", "referenceData/AN_genes_for_polII.bed")
orgDb <- org.Anidulans.eg.db
txDb <- TxDb.Anidulans.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")

file_polIISamples <- paste(polII_dataPath, "/", "sample_polII.list", sep = "")

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>%
  dplyr::select(-score) %>%
  dplyr::mutate(length = end - start)


geneDesc <- select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

##################################################################################
## process individual polII data
polIISampleList <- readr::read_tsv(file = file_polIISamples, col_names = c("id"),  comment = "#")

polII_info <- get_sample_information(exptInfoFile = file_exptInfo,
                                     samples = polIISampleList$id,
                                     dataPath = polII_dataPath,
                                     matrixSource = "normalizedmatrix")


polIICols <- list(
  exp = structure(polII_info$sampleId, names = polII_info$sampleId),
  is_expressed = structure(paste("is_expressed", ".", polII_info$sampleId, sep = ""),
                           names = polII_info$sampleId)
)



i <- 1

for (i in 1:nrow(polII_info)) {
  ## make gene level polII signal file for each sample
  polIIDf <- chipmine::preProcess_polII_expression(
    expMat = polII_info$polIIExpMat[i],
    sampleId = polII_info$sampleId[i],
    expFraction = 10,
    polIIExpFile = polII_info$polIIExpFile[i])
  
  # ## create 2kb - 2kb - 1kb matrix
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = polII_info$bwFile[i],
  #                                          regions = file_genes,
  #                                          signalName = polII_info$sampleId[i],
  #                                          extend = c(2000, 1000), w = 10,
  #                                          storeLocal = TRUE,
  #                                          localPath = polII_info$matFile[i],
  #                                          target_ratio = 0.4)
  
}


##################################################################################
## create polII data matrix

polIIMat <- get_polII_expressions(genesDf = geneSet, exptInfo = polII_info) %>% 
  # dplyr::select(-starts_with("is_expressed")) %>%
  dplyr::select(chr, start, end, geneId, strand, length, DESCRIPTION, everything())


readr::write_tsv(x = polIIMat, path = paste(polII_dataPath, "/polII_signal_matrix.tab", sep = ""))



##################################################################################
## polII signal percentile matrix

quantPoints <- c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1)


tmp <- quantile(x = polIIMat$AN0148_sCopy_OE_16h_polII_ChIPMix55_1, quantPoints, na.rm = T)


polIIQuantiles <- purrr::map_dfr(
  .x = polIICols$exp,
  .f = function(x){
    quantSummary <- as.list(round(quantile(x = polIIMat[[x]], quantPoints, na.rm = T), digits = 3))
    # quantSummary$sum <- sum(polIIMat[[x]])
    return(quantSummary)
  },
  .id = "sampleId"
)


readr::write_tsv(x = polIIQuantiles,
                 path = paste(polII_dataPath, "/polII_signal_quantiles.tab", sep = ""))


##################################################################################










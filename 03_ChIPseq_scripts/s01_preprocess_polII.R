library(chipmine)
library(org.AFumigatus.Af293.eg.db)
library(TxDb.Afumigatus.Af293.AspGD.GFF)
library(foreach)
library(doParallel)
library(here)
library(GGally)


rm(list = ls())
# cl <- makeCluster(4) #not to overload your computer
# registerDoParallel(cl)

##################################################################################

file_exptInfo <- here::here("data", "reference_data/sample_info.txt")

file_genes <- here::here("data", "reference_data/A_fumigatus_Af293_version_s03-m05-r12_genes_coding.bed")
orgDb <- org.AFumigatus.Af293.eg.db
txDb <- TxDb.Afumigatus.Af293.AspGD.GFF

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")

file_polIISamples <- paste(polII_dataPath, "/", "sample_polII.list", sep = "")

file_deeptolsMat <- paste(polII_dataPath, "/", "raw_count.deeptools.tab", sep = "")

##################################################################################

geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>%
  dplyr::select(-score) %>%
  dplyr::mutate(length = end - start)


geneDesc <- select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

polIISampleList <- readr::read_tsv(file = file_polIISamples, col_names = c("id"),  comment = "#")

polII_info <- get_sample_information(exptInfoFile = file_exptInfo,
                                     samples = polIISampleList$id,
                                     dataPath = polII_dataPath)


polIICols <- list(
  exp = structure(polII_info$sampleId, names = polII_info$sampleId),
  is_expressed = structure(paste("is_expressed", ".", polII_info$sampleId, sep = ""),
                           names = polII_info$sampleId)
)


##################################################################################
## process individual polII data

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
## process deeptools raw count matrix

deeptoolsMat <- suppressMessages(readr::read_tsv(file = file_deeptolsMat))

colnames(deeptoolsMat) <- stringr::str_replace(
  string = colnames(deeptoolsMat), pattern = "_bt2", replacement = "")

if(!setequal(
  paste(geneSet$chr, geneSet$start, geneSet$end, sep = ":"),
  paste(deeptoolsMat$chr, deeptoolsMat$start, deeptoolsMat$end, sep = ":"))){
  
  stop("geneSet genes and deeptools matrix genes does not match")
  
} else{
  rawCountMat <- dplyr::left_join(
    x = deeptoolsMat, y = geneSet, by = c("chr" , "start", "end")
  ) %>% 
    dplyr::select(geneId, polII_info$sampleId)
  
  readr::write_tsv(x = rawCountMat,
                   path = paste(polII_dataPath, "/polII_raw_counts.tab", sep = ""))
}















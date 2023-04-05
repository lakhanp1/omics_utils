library(chipmine)
library(here)
library(TxDb.Hsapiens.GRCh38p12.gencodev30.basic)

rm(list = ls())

file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")

TF_dataPath <- here::here("data", "ChIPseq_data")
txDb <- TxDb.Hsapiens.GRCh38p12.gencodev30.basic

##################################################################################
## get the sample details
tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")

geneSet <- GenomicFeatures::genes(x = txDb, columns = c("gene_id", "tx_id", "tx_name")) %>% 
  as.data.frame() %>% 
  dplyr::rename(
    chr = seqnames, geneId = gene_id
  )

i <- 6

for(i in 1:nrow(tfInfo)){
  
  peakAn <- narrowPeak_annotate(peakFile = tfInfo$narrowpeakFile[i],
                                txdb = txDb,
                                includeFractionCut = 0.7,
                                bindingInGene = FALSE,
                                promoterLength = 3000,
                                insideSkewToEndCut = 0.7,
                                reportPseudo = FALSE,
                                output = tfInfo$narrowpeakAnno[i])
  
  
  tfDf <- gene_level_peak_annotation(sampleId = tfInfo$sampleId[i],
                                     peakAnnotation = tfInfo$narrowpeakAnno[i],
                                     genesDf = geneSet,
                                     peakFile = tfInfo$narrowpeakFile[i],
                                     bwFile = tfInfo$bwFile[i],
                                     outFile = tfInfo$peakTargetFile[i],
                                     bindingInGene = FALSE)
  
  
}




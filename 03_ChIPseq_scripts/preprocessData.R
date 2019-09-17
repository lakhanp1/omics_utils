library(dplyr)
library(data.table)
library(chipmine)

rm(list = ls())


path = "E:/Chris_UM/Analysis/27_Brazil_data/"
setwd(path)

##################################################################################

exptInfoFile ="sampleInfo.txt"
TF_dataPath = "data"
polII_dataPath = "data"
geneCdsFile = "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r06_CDS_Unique.bed"
cdsUpRegionFile = ""


##################################################################################

polIIsampleFile = paste(polII_dataPath, "/", "polII_sample.list", sep = "")

polIISampleList = fread(file = polIIsampleFile, sep = "\t", header = F,
                        stringsAsFactors = F, col.names = c("id"), data.table = F)

polIISampleList = data.frame(id = c("An_untagged_20h_polII_1", "An_untagged_48h_polII_1"),
                             stringsAsFactors = F)

polII_info = get_sample_information(exptInfoFile = exptInfoFile,
                                    samples = polIISampleList$id,
                                    TF_path = TF_dataPath,
                                    polII_path = polII_dataPath)

## process all polII expression matrix
for(i in 1:nrow(polII_info)){
  
  polIIDf = preProcess_polII_expression(expMat = polII_info$polIIExpMat[i],
                                        title = polII_info$Sample_ID[i],
                                        expFraction = 10,
                                        polIIExpFile = polII_info$polIIExpFile[i])
  
  
}


##################################################################################
## process all TF macs2 results to merge everything

tfSampleFile = paste(TF_dataPath, "/", "TF_samples.list", sep = "")

tfSampleList = fread(file = tfSampleFile, sep = "\t", header = F,
                     stringsAsFactors = F, col.names = c("id"), data.table = F)

# tfSampleList = data.frame(id = c("An_laeA_20h_HA_1", "An_laeA_48h_HA_1"),
#                           stringsAsFactors = F)


tfInfo = get_sample_information(exptInfoFile = exptInfoFile,
                                samples = tfSampleList$id,
                                TF_path = TF_dataPath,
                                polII_path = polII_dataPath)

# i=1

for (i in 1:nrow(tfInfo)) {
  tfDf = preProcess_macs2_results(title = tfInfo$Sample_ID[i],
                                  peakAnnoFile = tfInfo$narrowpeakAnno[i],
                                  cdsFile = geneCdsFile,
                                  outFile = tfInfo$tfPeakFile[i],
                                  bindingInGene = FALSE)
}


## build peak region BED file for calculating peak region FPKM
## if gene has peak, use the peak region else use 1000bp upstream region
for (i in 1:nrow(tfInfo)) {
  genwise_peak_regions_bed(title = tfInfo$Sample_ID[i],
                           tfPeakFile = tfInfo$tfPeakFile[i],
                           cdsUpstreamFile = cdsUpRegionFile,
                           outFile = tfInfo$tfRegionsBed[i])
}

##################################################################################


## specific processing for samples where binding is seen in gene body
# tfSampleList = c("An_laeA_20h_HA", "An_laeA_48h_HA", "An_kdmB_20h_HA", "An_kdmB_48h_HA")
samplesWithBindingInGene = c("An_cclA_20h_HA_1", "An_cclA_48h_HA_1", "An_cclA_kdmA_del_20h_HA_1", "An_cclA_kdmA_del_48h_HA_1", "An_kdmA_20h_HA_1", "An_kdmA_48h_HA_1")


tfInfo = get_sample_information(exptInfoFile = exptInfoFile,
                                samples = samplesWithBindingInGene,
                                TF_path = TF_dataPath,
                                polII_path = polII_dataPath)


for (i in 1:nrow(tfInfo)) {
  tfDf = preProcess_macs2_results(title = tfInfo$Sample_ID[i],
                                  macs2Out = tfInfo$narrowpeakFile[i],
                                  peakAnnoFile = tfInfo$narrowpeakAnno[i],
                                  cdsFile = geneCdsFile,
                                  outFile = tfInfo$tfPeakFile[i],
                                  bindingInGene = TRUE)
}


## build peak region BED file for calculating peak region FPKM
for (i in 1:nrow(tfInfo)) {
  genwise_peak_regions_bed(title = tfInfo$Sample_ID[i],
                           tfPeakFile = tfInfo$tfPeakFile[i],
                           cdsUpstreamFile = cdsUpRegionFile,
                           outFile = tfInfo$tfRegionsBed[i])
}






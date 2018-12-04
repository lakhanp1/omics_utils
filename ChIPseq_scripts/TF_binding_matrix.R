library(dplyr)
library(data.table)
library(tibble)
library(chipmine)
library(GenomicRanges)
require(rtracklayer)
library(BSgenome.Afumigatus.AspGD.Af293)
library(ggplot2)


## This script generate a matrix of all TF 

rm(list = ls())

path = "E:/Chris_UM/Analysis/27_Brazil_data"
setwd(path)


outDir = "E:/Chris_UM/Analysis/27_Brazil_data/combined_analysis"


if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


exptInfoFile ="sampleInfo.txt"
TF_dataPath = "data"
file_genes = "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r06_CDS_Unique.bed"


sampleList = c("CREEHA_CONTROL2", "CREEHA_CONTROL3", "CREEHA_10MMAA2", "CREEHA_10MMAA3")


##################################################################################
## generate TF binding matrix based on the target genes

sampleList = c("CREEHA_CONTROL2", "CREEHA_10MMAA3")

genesDf = data.table::fread(file = file_genes, sep = "\t", header = F, select = c(4), col.names = c("gene"), stringsAsFactors = F)

## get the sample details
exptData = get_sample_information(exptInfoFile = exptInfoFile,
                                  samples = sampleList,
                                  TF_path = TF_dataPath,
                                  polII_path = TF_dataPath)


tfData = get_TF_binding_data(genesDf = genesDf, exptInfo = exptData) %>%
  dplyr::select(gene, starts_with("hasPeak"), starts_with("peakType"), starts_with("peakId"), starts_with("pval"), starts_with("enrichment"), starts_with("peakDist")) %>%
  dplyr::filter_at(.vars = vars(starts_with("peakType"), starts_with("tesPeakType")),
                   .vars_predicate = any_vars(!is.na(.)))


outFile = paste(outDir, "/peak_target_matrix_good_samples.tab", sep = "")

fwrite(x = tfData, file = outFile, sep = "\t", col.names = T, quote = F, na = "NA")


##################################################################################
## merge the peaks and generate the combination matrix for presence of peak

sampleList = c("CREEHA_CONTROL2", "CREEHA_10MMAA3")

## get the sample details
exptData = get_sample_information(exptInfoFile = exptInfoFile,
                                  samples = sampleList,
                                  TF_path = TF_dataPath,
                                  polII_path = TF_dataPath)

## generate the combinatorial binding matrix
mat = combinatorial_binding_matrix(sampleInfo = exptData,
                                   genome = BSgenome.Afumigatus.AspGD.Af293,
                                   summitSeqLen = 200)

sortCols = grep(pattern = "^pval", x = names(mat), value = T, perl = T)


newMat = dplyr::group_by_at(.tbl = mat, .vars = vars(starts_with("overlap"))) %>% 
  dplyr::arrange(desc(`pval(CREEHA_CONTROL2)`), desc(`pval(CREEHA_10MMAA3)`), .by_group = TRUE) %>% 
  dplyr::ungroup()


outFile = paste(outDir, "/peak_combinatorial_mat_201bp_seq.tab", sep = "")

data.table::fwrite(x = newMat, file = outFile, sep = "\t", col.names = T, quote = F, na = "NA")


# "pval(CREEHA_CONTROL2)", "pval(CREEHA_10MMAA3)" 
ggplot2::ggplot(data = mat, mapping = aes(x = `pval(CREEHA_10MMAA3)`)) +
  geom_histogram(binwidth = 2) +
  coord_cartesian(xlim = c(0, 80))











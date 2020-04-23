library(dplyr)
library(tidyr)
library(ChIPComp)


path = "E:/Chris_UM/Analysis/3_cAlbicanceVirulance/TBF1_analysis/TF_ChIP_analysis"
setwd(path)


compare = c("1254_MYC", "1255_MYC")

file_sampleInfo = "sample_info_tbf1.csv"


sampleInfo = fread(file = file_sampleInfo, sep = ",", header = T, stringsAsFactors = F) %>% 
  dplyr::select(sampleID = SampleID,
                condition = Condition,
                factor = Factor,
                ipReads = bamReads,
                ctReads = bamControl,
                peaks = bedPeaks)

conf = subset(sampleInfo, condition %in% compare)

design = as.data.frame(model.matrix(~condition, conf))

names(design) = c("(Intercept)", "condition")


countSet = makeCountSet(conf = conf, design = design, filetype="bam", species="other", binsize=10)

plot(countSet)

countSet = ChIPComp(countSet)

write.table(x = countSet$db, "clipboard", sep = "\t", col.names = T, row.names = F, quote = F)



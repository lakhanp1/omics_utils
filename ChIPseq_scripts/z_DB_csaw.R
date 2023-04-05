library(csaw)
library(dplyr)
library(tidyr)
library(DiffBind)
library(edgeR)
library(statmod)

rm(list = ls())


path = "E:/Chris_UM/Analysis/3_cAlbicanceVirulance/TBF1_analysis/TF_ChIP_analysis"
setwd(path)

file_sampleInfo = "sample_info_tbf1.csv"

sampleInfo = fread(file = file_sampleInfo, sep = ",", header = T, stringsAsFactors = F)

## DiffBind configuration
configDf = data.frame(RunParallel = TRUE,
                      # AnalysisMethod = DBA_DESEQ2,
                      ReportInit = "Tbf1_ChIP",
                      fragmentSize = 200)

##########################################################################
## load data from BAM files
fragLen = 200
winWid = 50
spaceWid = 40

param = readParam(minq = 30,
                  BPPARAM = SnowParam(workers = 6))

data = windowCounts(bam.files = sampleInfo$bamReads,
                    ext = fragLen,
                    spacing = spaceWid,
                    width = winWid,
                    param = param)

head(assay(data))
rowRanges(data)

##########################################################################
## Filter uninterested results
## gene body or global enrichment: polII ChIP data
## local enrichment or peak regions: TF ChIP data

## get the consensus peaks using DiffBind package
dataDiffBind = DiffBind::dba(sampleSheet = file_sampleInfo, config = configDf)
regions = dba.peakset(DBA = dataDiffBind, bRetrieve = T, DataType = DBA_DATA_GRANGES)

suppressWarnings(keep <- overlapsAny(rowRanges(data), regions))

filteredData = data[keep, ]


##########################################################################
## calculate normalization factors
## composition bias
binned = windowCounts(sampleInfo$bamReads, bin=TRUE, width=10000, param=param)

filteredData = normOffsets(binned, se.out = filteredData)

par(mfrow=c(3, 3), mar=c(5, 4, 2, 1.5))
adj.counts = cpm(asDGEList(binned), log=TRUE)

normfacs = filteredData$norm.factors

for (i in seq_len(length(sampleInfo$bamReads)-1)) {
  cur.x <- adj.counts[,1]
  cur.y <- adj.counts[,1+i]
  
  smoothScatter(x = (cur.x+cur.y)/2+6*log2(10),
                y = cur.x-cur.y,
                xlab = "A", ylab = "M", main = paste("1 vs", i+1)
                )
  
  all.dist = diff(log2(normfacs[c(i+1, 1)]))
  abline(h=all.dist, col="red")
  
}



## Efficiency bias
# keep2 = filterWindows(data = data, background = binned, type = "global")$filter > log2(5)
# 
# filtered.ac <- data[keep2,]
# normOffsets(filtered.ac, se.out = FALSE)
# 
# normOffsets(filteredData, se.out = FALSE)
# 
# filteredData$norm.factors = normOffsets(filtered.ac, se.out = FALSE)
# normfacs = filteredData$norm.factors

## trended bias
ac.demo2 = windowCounts(sampleInfo$bamReads, width=2000L, param=param)

filteredData = normOffsets(filteredData, type="loess", se.out = TRUE)

##########################################################################
## idnetify DB windows
y = asDGEList(object = filteredData, group = sampleInfo$Condition)

design = model.matrix(~ 0 + Condition, sampleInfo)
colnames(design) = levels(y$samples$group)

y = estimateDisp(y = y, design = design)
summary(y$trended.dispersion)

fit = glmQLFit(y, design, robust=TRUE)
summary(fit$var.post)


par(mfrow=c(1,2))
o = order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit)

results = glmQLFTest(fit, contrast=c(0, 0, -1, 1))
head(results$table)

rowData(filteredData) = cbind(rowData(filteredData), results$table)

##########################################################################
## correcting for multiple testing
olap = findOverlaps(regions, rowRanges(filteredData))

tabbroad = combineOverlaps(olap, results$table)


merged = mergeWindows(rowRanges(filteredData), tol=1000L)
merged$region

tabcom = combineTests(merged$id, results$table)

write.table(x = tabcom, "clipboard", sep = "\t", col.names = T, row.names = F, quote = F)

##########################################################################
## generate tabular output


##########################################################################
## data visualization





















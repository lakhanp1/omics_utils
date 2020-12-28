library(chipmine)
library(org.AFumigatus.Af293.eg.db)
library(cowplot)

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/topGO_functions.R")

outDir <- here::here("analysis", "02_ChIPseq_analysis", "03_diffbind", "topGO")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


outPrefix <- paste(outDir, "/topGO", sep = "")

##################################################################################
diffbindCompare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

file_goMap <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/annotation_resources/geneid2go.AFumigatus_Af293.topGO.map"
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")

TF_dataPath <- here::here("data", "TF_data")

sampleList <- c("CREEHA_CONTROL4", "CREEHA_CONTROL5", "CREEHA_10MMAA4", "CREEHA_10MMAA5")
file_diffbindTargets <- here::here("analysis", "02_ChIPseq_analysis",
                                   "01_peak_targets", "diffbind_allPeak_targets.tab")
orgDb <- org.AFumigatus.Af293.eg.db

##################################################################################

grp1 <- diffbindCompare[1]
grp2 <- diffbindCompare[2]
grp1Enrich = paste(grp1, ":enriched", sep = "")
grp2Enrich = paste(grp2, ":enriched", sep = "")
grp1Specific = paste(grp1, ":specific", sep = "")
grp2Specific = paste(grp2, ":specific", sep = "")


## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleList,
                                   dataPath = TF_dataPath,
                                   profileMatrixSuffix = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered"),
  FUN = function(x){ structure(paste(x, ".", sampleList, sep = ""), names = sampleList) },
  simplify = F, USE.NAMES = T
)

diffbindTargets <- suppressMessages(readr::read_tsv(file = file_diffbindTargets, col_names = T))
diffbindTargets <- dplyr::filter(diffbindTargets, pvalGood.all > 0)

dplyr::group_by(diffbindTargets, categoryDiffbind) %>% 
  dplyr::summarise(n = n())

nodeSize <- 5

##################################################################################

## peaks with stronger binding in CREEHA_CONTROL sample
outTf1Enrich <- paste(outPrefix, ".", grp1, "_enriched", sep = "")

tf1Enrich <- dplyr::filter(.data = diffbindTargets,
                                  categoryDiffbind %in% c(grp1Enrich, grp1Specific))

nrow(tf1Enrich)

tf1EnrichGo <- topGO_enrichment(genes = unique(tf1Enrich$geneId),
                          goMapFile = file_goMap, goNodeSize = nodeSize)

tf1EnrichGo$condition <- grp1Enrich

readr::write_tsv(x = tf1EnrichGo, path = paste(outTf1Enrich, ".tab", sep = ""))

plotTitle <- paste(grp1, ": increased binding signal (n = ", length(unique(tf1Enrich$geneId)), ")",
                   sep = "")

pt_tf1EnrichGo <- enrichment_scatter(df = tf1EnrichGo, title = plotTitle)

##################################################################################
## peaks specific to CREEHA_CONTROL sample
outTf1Specific <- paste(outPrefix, ".", grp1, "_specific", sep = "")

tf1Specific <- dplyr::filter(.data = diffbindTargets,
                             categoryDiffbind %in% c(grp1Specific))

nrow(tf1Specific)

tf1SpecificGo <- topGO_enrichment(genes = unique(tf1Specific$geneId),
                          goMapFile = file_goMap, goNodeSize = nodeSize)

tf1SpecificGo$condition <- grp1Specific

readr::write_tsv(x = tf1SpecificGo, path = paste(outTf1Specific, ".tab", sep = ""))

plotTitle <- paste(grp1, ": specific peaks (n = ", length(unique(tf1Specific$geneId)), ")",
                   sep = "")

pt_tf1SpecificGo <- enrichment_scatter(df = tf1SpecificGo, title = plotTitle)

##################################################################################
## peaks with stronger binding in CREEHA_10MMAA3 sample
outTf2Enrich <- paste(outPrefix, ".", grp2, "_enriched", sep = "")

tf2Enrich <- dplyr::filter(.data = diffbindTargets,
                                  categoryDiffbind %in% c(grp2Enrich, grp2Specific))
nrow(tf2Enrich)

tf2EnrichGo <- topGO_enrichment(genes = unique(tf2Enrich$geneId),
                          goMapFile = file_goMap, goNodeSize = nodeSize)

tf2EnrichGo$condition <- grp2Enrich

readr::write_tsv(x = tf2EnrichGo, path = paste(outTf2Enrich, ".tab", sep = ""))

plotTitle <- paste(grp2, ": increased binding signal (n = ", length(unique(tf2Enrich$geneId)), ")",
                   sep = "")

pt_tf2EnrichGo <- enrichment_scatter(df = tf2EnrichGo, title = plotTitle)

##################################################################################
## peaks with stronger binding in CREEHA_10MMAA3 sample

outTf2Specific <- paste(outPrefix, ".", grp2, "_specific", sep = "")

tf2Specific <- dplyr::filter(.data = diffbindTargets,
                             categoryDiffbind %in% c(grp2Specific))
nrow(tf2Specific)

tf2SpecificGo <- topGO_enrichment(genes = unique(tf2Specific$geneId),
                          goMapFile = file_goMap, goNodeSize = nodeSize)

tf2SpecificGo$condition <- grp2Specific

readr::write_tsv(x = tf2SpecificGo, path = paste(outTf2Specific, ".tab", sep = ""))

plotTitle <- paste(grp2, ": specific peaks (n = ", length(unique(tf2Specific$geneId)), ")",
                   sep = "")

pt_tf2SpecificGo <- enrichment_scatter(df = tf2SpecificGo, title = plotTitle)

##################################################################################
## common targets for CREEHA_CONTROL2 and CREEHA_10MMAA3 samples with no DiffBind change

outcommonNoDiff <- paste(outPrefix, ".", "common_noDiff", sep = "")

commonNoDiff <- dplyr::filter(.data = diffbindTargets, categoryDiffbind == "common")

nrow(commonNoDiff)

commonNoDiffGo <- topGO_enrichment(genes = unique(commonNoDiff$geneId),
                             goMapFile = file_goMap, goNodeSize = nodeSize)

commonNoDiffGo$condition <- "common_noDiff"

readr::write_tsv(x = commonNoDiffGo, path = paste(outcommonNoDiff, ".tab", sep = ""))

plotTitle <- paste("common peaks between ", grp1, " and ",
                   grp2, "\nwith similar signal (n = ", length(unique(commonNoDiff$geneId)), ")",
                   sep = "")
pt_commonNoDiffGo <- enrichment_scatter(df = commonNoDiffGo, title = plotTitle)


##################################################################################
## common targets for CREEHA_CONTROL2 and CREEHA_10MMAA3 samples with no DiffBind change

outCommonBound <- paste(outPrefix, ".", "common_bound", sep = "")

commonBound <- dplyr::filter(.data = diffbindTargets,
                             categoryDiffbind %in% c(grp1Enrich, "common", grp2Enrich))

nrow(commonBound)

commonBoundGo <- topGO_enrichment(genes = unique(commonBound$geneId),
                                  goMapFile = file_goMap, goNodeSize = nodeSize)

commonBoundGo$condition <- "common_bound"

readr::write_tsv(x = commonBoundGo, path = paste(outCommonBound, ".tab", sep = ""))

plotTitle <- paste("common peaks:", grp1, " and ",
                   grp2, " (n = ", length(unique(commonBound$geneId)), ")",
                   sep = "")
pt_commonBoundGo <- enrichment_scatter(df = commonBoundGo, title = plotTitle)

##################################################################################

aligned_plots <- align_plots(
  plotlist = list(
    tf1Enrich = pt_tf1EnrichGo,
    tf1Specific = pt_tf1SpecificGo,
    tf2Enrich = pt_tf2EnrichGo,
    tf2Specific = pt_tf2SpecificGo,
    commonNoDiff = pt_commonNoDiffGo,
    commonBound = pt_commonBoundGo
  ),
  align = "v")

pdf(file = paste(outPrefix, ".scatterplot.pdf", sep = ""),
    width = 15, height = 14, onefile = T, pointsize = 18)
ggdraw(aligned_plots$tf1Enrich)
ggdraw(aligned_plots$tf1Specific)
ggdraw(aligned_plots$tf2Enrich)
ggdraw(aligned_plots$tf2Specific)
ggdraw(aligned_plots$commonNoDiff)
ggdraw(aligned_plots$commonBound)
dev.off()


# wd <- (min(max(nchar(as.character(goData$Term))), 80) * 30) * 1.5
wd <- 3500
ht <- 4000
# res = max(min(wd, ht) / 12, 200)
res <- 300
# ht <- max(nrow(tf1EnrichGo) * 80, 1500)
png(filename = paste(outTf1Enrich, ".scatter.png", sep = ""),
    width = wd, height = ht, res = res)
ggdraw(aligned_plots$tf1Enrich)
dev.off()


# ht <- max(nrow(tf1SpecificGo) * 80, 1500)
png(filename = paste(outTf1Specific, ".scatter.png", sep = ""),
    width = wd, height = ht, res = res)
ggdraw(aligned_plots$tf1Specific)
dev.off()

# ht <- max(nrow(tf2EnrichGo) * 80, 1500)
png(filename = paste(outTf2Enrich, ".scatter.png", sep = ""),
    width = wd, height = ht, res = res)
ggdraw(aligned_plots$tf2Enrich)
dev.off()


# ht <- max(nrow(tf2SpecificGo) * 80, 1500)
png(filename = paste(outTf2Specific, ".scatter.png", sep = ""),
    width = wd, height = ht, res = res)
ggdraw(aligned_plots$tf2Specific)
dev.off()


# ht <- max(nrow(commonNoDiffGo) * 80, 1500)
png(filename = paste(outcommonNoDiff, ".scatter.png", sep = ""),
    width = wd, height = ht, res = res)
ggdraw(aligned_plots$commonNoDiff)
dev.off()

# ht <- max(nrow(commonBoundGo) * 80, 1500)
png(filename = paste(outCommonBound, ".scatter.png", sep = ""),
    width = wd, height = ht, res = res)
ggdraw(aligned_plots$commonBound)
dev.off()


##################################################################################





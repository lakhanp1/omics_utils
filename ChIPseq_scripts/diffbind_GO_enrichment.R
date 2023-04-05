suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(cowplot))


rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/s01_enrichment_functions.R")

outDir <- here::here("analysis", "11_aflR_AN7820_analysis", "a04_diffbind")

analysisName <- "AflR_diffbind"
outPrefix <- paste(outDir, "/", analysisName, ".topGO", sep = "")

## the denominator or WT in log2(fold_change) should be second
col_compare <- "diffbind_condition"
diffbindCompare <- c("high_Xylose", "low_Xylose")

##################################################################################

file_topGO <- "E:/Chris_UM/Database/A_Nidulans/annotation_resources/geneid2go.ANidulans.topGO.map"

file_diffbindInfo <- paste(outDir, "/diffbind_info.txt", sep = "")
file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")
file_diffbindRes <- paste(outDir, "/AflR_diffbind.annotation.filtered.tab", sep = "")

TF_dataPath <- here::here("data", "TF_data")


orgDb <- org.Anidulans.FGSCA4.eg.db
keggOrg <- 'ani'									                  ## KEGG organism code
col_degIdKeytype <- "GID"               ## org.db keytype of DEG file geneId
col_kegg <- "NCBI_ID"																## org.db column for NCBI ID
col_gsea <- "GID"																## org.db column to use for gsea analysis
col_topGO <- "GID"															## org.db keytype of topGO map file geneId
col_geneName <- "GENE_NAME"													## org.db column for Gene name

##################################################################################

if(!dir.exists(outDir)){dir.create(path = outDir)}

grp1 <- diffbindCompare[1]
grp2 <- diffbindCompare[2]
grp1Enrich = paste(grp1, ":enriched", sep = "")
grp2Enrich = paste(grp2, ":enriched", sep = "")
grp1Specific = paste(grp1, ":specific", sep = "")
grp2Specific = paste(grp2, ":specific", sep = "")

diffbindInfo <- suppressMessages(readr::read_tsv(file = file_diffbindInfo))

## get the sample details
exptData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = diffbindInfo$sampleId,
  dataPath = TF_dataPath
)

exptData <- dplyr::left_join(x = exptData, y = diffbindInfo, by = "sampleId") %>% 
  dplyr::mutate(
    sampleName = sampleLabel
  )

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))


diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffbindRes)) %>% 
  dplyr::filter(pvalGood.all != 0)

# glimpse(diffbindRes)

dplyr::group_by(diffbindRes, categoryDiffbind) %>% 
  dplyr::summarise(n = n())

nodeSize <- 5
##################################################################################

## peaks with stronger binding in grp1
outDiffbindUp <- paste(outPrefix, ".up", sep = "")

diffbindUp <- dplyr::filter(diffbindRes, diffBind == "up")

nrow(diffbindUp)

## topGO GO enrichment
diffbindUp_Go <- topGO_enrichment(
  goMapFile = file_topGO,
  genes = as.character(na.exclude(unique(diffbindUp$geneId))),
  type = "BP", goNodeSize = 5,
  orgdb = orgDb, keytype = col_degIdKeytype,
  topgoColumn = col_topGO, geneNameColumn = col_geneName
)


diffbindUp_Go$condition <- "up"

readr::write_tsv(x = diffbindUp_Go, file = paste(outDiffbindUp, ".tab", sep = ""))

plotTitle <- paste(grp1, ": increased binding signal (n = ", diffbindUp_Go$inputSize[1], ")",
                   sep = "")

pt_diffbindUp_Go <- enrichment_scatter(df = diffbindUp_Go, title = plotTitle)

##################################################################################

## peaks with stronger binding in grp2
outDiffbindDown <- paste(outPrefix, ".down", sep = "")

diffbindDown <- dplyr::filter(diffbindRes, diffBind == "down")

nrow(diffbindDown)

## topGO GO enrichment
diffbindDown_Go <- topGO_enrichment(
  goMapFile = file_topGO,
  genes = as.character(na.exclude(unique(diffbindDown$geneId))),
  type = "BP", goNodeSize = 5,
  orgdb = orgDb, keytype = col_degIdKeytype,
  topgoColumn = col_topGO, geneNameColumn = col_geneName
)


diffbindDown_Go$condition <- "down"

readr::write_tsv(x = diffbindDown_Go, file = paste(outDiffbindDown, ".tab", sep = ""))

plotTitle <- paste(grp2, ": increased binding signal (n = ", diffbindDown_Go$inputSize[1], ")",
                   sep = "")

pt_diffbindDown_Go <- enrichment_scatter(df = diffbindDown_Go, title = plotTitle)

##################################################################################
# ## peaks with stronger binding in CREEHA_CONTROL sample
# outTf1Enrich <- paste(outPrefix, ".", grp1, "_enriched", sep = "")
# 
# tf1Enrich <- dplyr::filter(.data = diffbindRes,
#                                   categoryDiffbind %in% c(grp1Enrich, grp1Specific))
# 
# nrow(tf1Enrich)
# 
# tf1EnrichGo <- topGO_enrichment(genes = unique(tf1Enrich$geneId),
#                           goMapFile = file_topGO, goNodeSize = nodeSize)
# 
# tf1EnrichGo$condition <- grp1Enrich
# 
# readr::write_tsv(x = tf1EnrichGo, path = paste(outTf1Enrich, ".tab", sep = ""))
# 
# plotTitle <- paste(grp1, ": increased binding signal (n = ", length(unique(tf1Enrich$geneId)), ")",
#                    sep = "")
# 
# pt_tf1EnrichGo <- enrichment_scatter(df = tf1EnrichGo, title = plotTitle)
# 
# ##################################################################################
# ## peaks specific to CREEHA_CONTROL sample
# outTf1Specific <- paste(outPrefix, ".", grp1, "_specific", sep = "")
# 
# tf1Specific <- dplyr::filter(.data = diffbindRes,
#                              categoryDiffbind %in% c(grp1Specific))
# 
# nrow(tf1Specific)
# 
# tf1SpecificGo <- topGO_enrichment(genes = unique(tf1Specific$geneId),
#                           goMapFile = file_topGO, goNodeSize = nodeSize)
# 
# tf1SpecificGo$condition <- grp1Specific
# 
# readr::write_tsv(x = tf1SpecificGo, path = paste(outTf1Specific, ".tab", sep = ""))
# 
# plotTitle <- paste(grp1, ": specific peaks (n = ", length(unique(tf1Specific$geneId)), ")",
#                    sep = "")
# 
# pt_tf1SpecificGo <- enrichment_scatter(df = tf1SpecificGo, title = plotTitle)
# 
# ##################################################################################
# ## peaks with stronger binding in CREEHA_10MMAA3 sample
# outTf2Enrich <- paste(outPrefix, ".", grp2, "_enriched", sep = "")
# 
# tf2Enrich <- dplyr::filter(.data = diffbindRes,
#                                   categoryDiffbind %in% c(grp2Enrich, grp2Specific))
# nrow(tf2Enrich)
# 
# tf2EnrichGo <- topGO_enrichment(genes = unique(tf2Enrich$geneId),
#                           goMapFile = file_topGO, goNodeSize = nodeSize)
# 
# tf2EnrichGo$condition <- grp2Enrich
# 
# readr::write_tsv(x = tf2EnrichGo, path = paste(outTf2Enrich, ".tab", sep = ""))
# 
# plotTitle <- paste(grp2, ": increased binding signal (n = ", length(unique(tf2Enrich$geneId)), ")",
#                    sep = "")
# 
# pt_tf2EnrichGo <- enrichment_scatter(df = tf2EnrichGo, title = plotTitle)
# 
# ##################################################################################
# ## peaks with stronger binding in CREEHA_10MMAA3 sample
# 
# outTf2Specific <- paste(outPrefix, ".", grp2, "_specific", sep = "")
# 
# tf2Specific <- dplyr::filter(.data = diffbindRes,
#                              categoryDiffbind %in% c(grp2Specific))
# nrow(tf2Specific)
# 
# tf2SpecificGo <- topGO_enrichment(genes = unique(tf2Specific$geneId),
#                           goMapFile = file_topGO, goNodeSize = nodeSize)
# 
# tf2SpecificGo$condition <- grp2Specific
# 
# readr::write_tsv(x = tf2SpecificGo, path = paste(outTf2Specific, ".tab", sep = ""))
# 
# plotTitle <- paste(grp2, ": specific peaks (n = ", length(unique(tf2Specific$geneId)), ")",
#                    sep = "")
# 
# pt_tf2SpecificGo <- enrichment_scatter(df = tf2SpecificGo, title = plotTitle)
# 
# ##################################################################################
# ## common targets for CREEHA_CONTROL2 and CREEHA_10MMAA3 samples with no DiffBind change
# 
# outcommonNoDiff <- paste(outPrefix, ".", "common_noDiff", sep = "")
# 
# commonNoDiff <- dplyr::filter(.data = diffbindRes, categoryDiffbind == "common")
# 
# nrow(commonNoDiff)
# 
# commonNoDiffGo <- topGO_enrichment(genes = unique(commonNoDiff$geneId),
#                              goMapFile = file_topGO, goNodeSize = nodeSize)
# 
# commonNoDiffGo$condition <- "common_noDiff"
# 
# readr::write_tsv(x = commonNoDiffGo, path = paste(outcommonNoDiff, ".tab", sep = ""))
# 
# plotTitle <- paste("common peaks between ", grp1, " and ",
#                    grp2, "\nwith similar signal (n = ", length(unique(commonNoDiff$geneId)), ")",
#                    sep = "")
# pt_commonNoDiffGo <- enrichment_scatter(df = commonNoDiffGo, title = plotTitle)
# 
# 
# ##################################################################################
# ## common targets for CREEHA_CONTROL2 and CREEHA_10MMAA3 samples with no DiffBind change
# 
# outCommonBound <- paste(outPrefix, ".", "common_bound", sep = "")
# 
# commonBound <- dplyr::filter(.data = diffbindRes,
#                              categoryDiffbind %in% c(grp1Enrich, "common", grp2Enrich))
# 
# nrow(commonBound)
# 
# commonBoundGo <- topGO_enrichment(genes = unique(commonBound$geneId),
#                                   goMapFile = file_topGO, goNodeSize = nodeSize)
# 
# commonBoundGo$condition <- "common_bound"
# 
# readr::write_tsv(x = commonBoundGo, path = paste(outCommonBound, ".tab", sep = ""))
# 
# plotTitle <- paste("common peaks:", grp1, " and ",
#                    grp2, " (n = ", length(unique(commonBound$geneId)), ")",
#                    sep = "")
# pt_commonBoundGo <- enrichment_scatter(df = commonBoundGo, title = plotTitle)
# 
##################################################################################

aligned_plots <- align_plots(
  plotlist = list(
    # tf1Enrich = pt_tf1EnrichGo,
    # tf1Specific = pt_tf1SpecificGo,
    # tf2Enrich = pt_tf2EnrichGo,
    # tf2Specific = pt_tf2SpecificGo,
    # commonNoDiff = pt_commonNoDiffGo,
    # commonBound = pt_commonBoundGo
    diffbindUp = pt_diffbindUp_Go,
    diffbindDown = pt_diffbindDown_Go
  ),
  align = "v")

pdf(file = paste(outPrefix, ".scatterplot.pdf", sep = ""),
    width = 10, height = 6, onefile = T, pointsize = 18)
# ggdraw(aligned_plots$tf1Enrich)
# ggdraw(aligned_plots$tf1Specific)
# ggdraw(aligned_plots$tf2Enrich)
# ggdraw(aligned_plots$tf2Specific)
# ggdraw(aligned_plots$commonNoDiff)
# ggdraw(aligned_plots$commonBound)
ggdraw(aligned_plots$diffbindUp)
ggdraw(aligned_plots$diffbindDown)
dev.off()


# # wd <- (min(max(nchar(as.character(goData$Term))), 80) * 30) * 1.5
# wd <- 3500
# ht <- 4000
# # res = max(min(wd, ht) / 12, 200)
# res <- 300
# # ht <- max(nrow(tf1EnrichGo) * 80, 1500)
# png(filename = paste(outTf1Enrich, ".scatter.png", sep = ""),
#     width = wd, height = ht, res = res)
# ggdraw(aligned_plots$tf1Enrich)
# dev.off()
# 
# 
# # ht <- max(nrow(tf1SpecificGo) * 80, 1500)
# png(filename = paste(outTf1Specific, ".scatter.png", sep = ""),
#     width = wd, height = ht, res = res)
# ggdraw(aligned_plots$tf1Specific)
# dev.off()
# 
# # ht <- max(nrow(tf2EnrichGo) * 80, 1500)
# png(filename = paste(outTf2Enrich, ".scatter.png", sep = ""),
#     width = wd, height = ht, res = res)
# ggdraw(aligned_plots$tf2Enrich)
# dev.off()
# 
# 
# # ht <- max(nrow(tf2SpecificGo) * 80, 1500)
# png(filename = paste(outTf2Specific, ".scatter.png", sep = ""),
#     width = wd, height = ht, res = res)
# ggdraw(aligned_plots$tf2Specific)
# dev.off()
# 
# 
# # ht <- max(nrow(commonNoDiffGo) * 80, 1500)
# png(filename = paste(outcommonNoDiff, ".scatter.png", sep = ""),
#     width = wd, height = ht, res = res)
# ggdraw(aligned_plots$commonNoDiff)
# dev.off()
# 
# # ht <- max(nrow(commonBoundGo) * 80, 1500)
# png(filename = paste(outCommonBound, ".scatter.png", sep = ""),
#     width = wd, height = ht, res = res)
# ggdraw(aligned_plots$commonBound)
# dev.off()


##################################################################################





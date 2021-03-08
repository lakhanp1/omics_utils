suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(DiffBind))


## this script:
## 1) run DiffBind to perform differential binding analysis
## 2) add peak data for individual sample
## 3) prepare combined report
## 4) add peak target gene information for each best peaksets

rm(list = ls())

outDir <- here::here("analysis", "11_aflR_AN7820_analysis", "a04_diffbind")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

## the denominator or WT in log2(fold_change) should be second
col_compare <- "Condition"
compare <- c("AN7820_sCopy_OE_minor_Xyl", "AN7820_sCopy_OE")

analysisName <- "AflR_diffbind"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

fdr_cut <- 0.05
lfc_cut <- 1
up_cut <- lfc_cut
down_cut <- lfc_cut * -1

##################################################################################
orgDb <- org.Anidulans.FGSCA4.eg.db

file_diffbindInfo <- paste(outDir, "/diffbind_info.txt", sep = "")

file_exptInfo <- here::here("data", "reference_data", "sample_info.txt")

TF_dataPath <- here::here("data", "TF_data")

##################################################################################
diffbindInfo <- suppressMessages(readr::read_tsv(file = file_diffbindInfo))

## get the sample details
exptData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = diffbindInfo$sampleId,
  dataPath = TF_dataPath
)

exptData$bam <- paste(here::here("data"), "/bams/", exptData$sampleId, "_bt2.bam", sep = "")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfCols <- sapply(
  c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
    "peakAnnotation", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
    "peakCoverage", "pvalFiltered", "summitSeq", "overlap", "targetGene"),
  FUN = function(x){ structure(paste(x, ".", exptData$sampleId, sep = ""), names = exptData$sampleId) },
  simplify = F, USE.NAMES = T)

diffbindInfo <- dplyr::select(
  exptData,
  SampleID = sampleId,
  Condition = condition,
  Replicate = rep,
  bamReads = bam,
  Peaks = peakFile,
  repRank
) %>% 
  dplyr::mutate(PeakCaller = "narrow")

if(is.null(diffbindInfo$pval_cutoff)){
  diffbindInfo$pval_cutoff <- 20
}

## ensure that reference level is same as compare[2] in the factor
diffbindInfo <- diffbindInfo %>% dplyr::mutate(
  !!col_compare := forcats::fct_relevel(.f = !!sym(col_compare), compare[2], compare[1])
) %>% 
  dplyr::arrange(!!sym(col_compare))


grp1 <- compare[1]        ## numerator
grp2 <- compare[2]        ## denominator
grp1Index <- which(diffbindInfo[[col_compare]] == grp1)
grp2Index <- which(diffbindInfo[[col_compare]] == grp2)
grp1Samples <- diffbindInfo$SampleID[grp1Index]
grp2Samples <- diffbindInfo$SampleID[grp2Index]
grp1SpecificOcc = paste(grp1, ":specific", sep = "")
grp2SpecificOcc = paste(grp2, ":specific", sep = "")

bestGrp1Id <- diffbindInfo$SampleID[diffbindInfo[[col_compare]] == grp1 & diffbindInfo$repRank == 1]
bestGrp2Id <- diffbindInfo$SampleID[diffbindInfo[[col_compare]] == grp2 & diffbindInfo$repRank == 1]

groupCols <- sapply(
  X = c("peakCall", "pvalGood"),
  FUN = function(x){ structure(paste(x, ".", compare, sep = ""), names = compare) },
  simplify = F, USE.NAMES = T
)

##################################################################################
## diffBind analysis
# 
# dataDba <- DiffBind::dba(sampleSheet = diffbindInfo)
# 
# plot(dataDba)
# 
# countsDba <- DiffBind::dba.count(DBA = dataDba)
# 
# plot(countsDba)
# 
# normDba <- dba.normalize(countsDba)
# 
# contrastDba <- DiffBind::dba.contrast(
#   DBA = normDba, design = "~Condition", minMembers = 2,
#   contrast = c(col_compare, compare)
# )
# 
# dba.show(contrastDba, bContrasts=TRUE)
# 
# diffDba <- DiffBind::dba.analyze(DBA = contrastDba, method = DBA_DESEQ2)
# 
# DiffBind::dba.save(
#   DBA = diffDba, file = paste(analysisName, ".dba", sep = ""),
#   dir = outDir, pre = ""
# )


## load DBA object
diffDba <- DiffBind::dba.load(file = paste(analysisName, ".dba", sep = ""), dir = outDir, pre = "")

dba.show(diffDba, bContrasts=TRUE)

plot(diffDba, contrast=1)

diffDf <- DiffBind::dba.report(
  DBA = diffDba, bFlip = TRUE, th = 1,
  bCalledDetail = T, DataType = DBA_DATA_FRAME
)

diffDf <- dplyr::mutate(diffDf, region = paste(Chr, ":", Start, "-", End, sep = ""))

# dba.plotMA(diffDba)
# dba.plotVolcano(diffDba)
# dba.plotBox(diffDba)
# dba.overlap(DBA = diffDba, mask = diffDba$masks$All, mode = DBA_OLAP_PEAKS)
# dba.peakset(DBA = diffDba, consensus = )
# htDt <- dba.plotHeatmap(diffDba, maxSites = diffDba$totalMerged,
#                         score = DBA_SCORE_TMM_MINUS_FULL, correlations=FALSE, th = 1)

##################################################################################
## create report table

## add individual peak information columns
diffGr <- DiffBind::dba.report(DBA = diffDba, th = 1)
diffGr <- sort(diffGr)
mcols(diffGr)$name <- paste("peak_region", 1:length(diffGr), sep = "_")

rtracklayer::export.bed(
  object = diffGr, format = "bed", 
  con = paste(outPrefix, ".merged_regions.bed", sep = "")
)

## add overlapping peaks from each sample
diffRes <- chipmine::combinatorial_binding_matrix(sampleInfo = exptData, peakRegions = diffGr)

## count of samples which showed macs2 peak under condition1
diffRes[[unname(groupCols$peakCall[grp1])]] <- purrr::pmap_int(
  .l = dplyr::select(diffRes, unname(tfCols$overlap[grp1Samples])),
  .f = sum, na.rm = TRUE
)

## count of samples which showed macs2 peak under condition2
diffRes[[unname(groupCols$peakCall[grp2])]] <- purrr::pmap_int(
  .l = dplyr::select(diffRes, unname(tfCols$overlap[grp2Samples])),
  .f = sum, na.rm = TRUE
)


dplyr::group_by_at(diffRes, .vars = vars(starts_with("peakCall"))) %>% 
  dplyr::summarise(n = n())

## add diffBind status and peak occupancy status
diffData <- diffRes %>% 
  dplyr::mutate(
    diffBind = dplyr::case_when(
      Fold >= up_cut & FDR <= fdr_cut ~ "up",
      Fold <= down_cut & FDR <= fdr_cut ~ "down",
      TRUE ~ "noDiff"
    ),
    peakOccupancy = dplyr::case_when(
      !!sym(groupCols$peakCall[grp1]) >= 2 & !!sym(groupCols$peakCall[grp2]) >= 2 ~ "common",
      !!sym(groupCols$peakCall[grp1]) >= 2 & !!sym(groupCols$peakCall[grp2]) < 2 ~ 
        !! grp1SpecificOcc,
      !!sym(groupCols$peakCall[grp1]) < 2 & !!sym(groupCols$peakCall[grp2]) >= 2 ~ 
        !! grp2SpecificOcc,
      TRUE ~ "no_consensus"
    )
  )

## add counts for number of samples showing macs2 pval > cutoff: group1
diffData[[ groupCols$pvalGood[grp1] ]] <- purrr::pmap_int(
  .l = dplyr::select(diffData, unname(tfCols$peakPval[grp1Samples])),
  .f = function(...){
    sum(c(...) >= diffbindInfo$pval_cutoff[grp1Index], na.rm = TRUE)
  }
)

## add counts for number of samples showing macs2 pval > cutoff: group2
diffData[[ groupCols$pvalGood[grp2] ]] <- purrr::pmap_int(
  .l = dplyr::select(diffData, unname(tfCols$peakPval[grp2Samples])),
  .f = function(...){
    sum(c(...) >= diffbindInfo$pval_cutoff[grp2Index], na.rm = TRUE)
  }
)

## add counts for total number of samples showing macs2 pval > cutoff
diffData$pvalGood.all <- purrr::pmap_int(
  .l = dplyr::select(diffData, unname(groupCols$pvalGood)),
  .f = sum, na.rm = TRUE
)


dplyr::group_by(diffData, diffBind, peakOccupancy) %>% 
  dplyr::summarise(n = n())

## assign peak category based on diffbind and peakOccupancy
diffData <- diffData %>% 
  dplyr::mutate(
    # categoryDiffbind = dplyr::case_when(
    #   diffBind == "up" | peakOccupancy == grp1SpecificOcc ~ paste(grp1, ":enriched", sep = ""),
    #   diffBind == "down" | peakOccupancy == grp2SpecificOcc ~ paste(grp2, ":enriched", sep = ""),
    #   diffBind == "noDiff" & peakOccupancy == "common" ~ "common",
    #   TRUE ~ "NA"
    # ),
    categoryDiffbind = dplyr::case_when(
      peakOccupancy == grp1SpecificOcc ~ grp1SpecificOcc,
      peakOccupancy == grp2SpecificOcc ~ grp2SpecificOcc,
      diffBind == "up" ~ paste(grp1, ":enriched", sep = ""),
      diffBind == "down" ~ paste(grp2, ":enriched", sep = ""),
      diffBind == "noDiff" & peakOccupancy == "common" ~ "common",
      TRUE ~ "NA"
    )
  )

diffData <- diffData %>% 
  dplyr::select(
    seqnames, start, end, name, starts_with("Conc"), Fold, p.value, FDR, diffBind,
    peakOccupancy, categoryDiffbind, starts_with("peakCall."), starts_with("pvalGood."),
    everything(), -starts_with("peakEnrichment."), -starts_with("overlap")
  )


readr::write_tsv(x = diffData, file = paste(outPrefix, ".all.tab", sep = ""))

##################################################################################
## diffbind report with target genes

# 
# ## prepare txIds excluding rRNA and tRNA transcripts
# geneInfo <- AnnotationDbi::select(x = orgDb,
#                                   keys = keys(orgDb, keytype = "GID"),
#                                   columns = c("GENE_NAME", "TYPE"),
#                                   keytype = "GID") %>% 
#   dplyr::rename(geneId = GID)
# 
# geneInfo %>% dplyr::filter(!grepl(pattern = "ORF\\|", x = TYPE, perl = TRUE)) %>% 
#   dplyr::select(geneId, GENE_NAME, TYPE) %>%
#   dplyr::count(TYPE)
# 
# geneInfo <- dplyr::filter(geneInfo, !grepl(pattern = "(rRNA|tRNA)\\|", x = TYPE, perl = TRUE))
# 
# txInfo <- AnnotationDbi::select(x = txDb, keys = geneInfo$geneId,
#                                 columns = "TXID", keytype = "GENEID")
# 
# ## annotate DiffBind regions
# diffGrAn <- annotate_ranges(peaks = diffGr, txdb = txDb, promoterLength = 500, txIds = txInfo$TXID)


## import peak annotation for sample1 and add new column with pval cutoff pass result
s1Targets <- import_peak_annotation(
  sampleId = bestGrp1Id,
  peakAnnoFile = exptDataList[[bestGrp1Id]]$peakAnno,
  columns = c("peakId", "geneId", "peakPosition", "peakAnnotation", "peakDist")
) %>% 
  dplyr::filter(!is.na(geneId)) %>% 
  dplyr::filter(!grepl(pattern = "pseudo_", x = !!sym(tfCols$peakAnnotation[bestGrp1Id]))) %>% 
  dplyr::rename(!! sym(tfCols$targetGene[bestGrp1Id]) := geneId)


## import peak annotation for sample2 and add new column with pval cutoff pass result
s2Targets <- import_peak_annotation(
  sampleId = bestGrp2Id,
  peakAnnoFile = exptDataList[[bestGrp2Id]]$peakAnno,
  columns = c("peakId", "geneId", "peakPosition", "peakAnnotation", "peakDist")
) %>% 
  dplyr::filter(!is.na(geneId)) %>% 
  dplyr::filter(!grepl(pattern = "pseudo_", x = !!sym(tfCols$peakAnnotation[bestGrp2Id]))) %>% 
  dplyr::rename(!! sym(tfCols$targetGene[bestGrp2Id]) := geneId)



## combine diffbind results with peak target annotation information
diffAnn <- diffData %>% 
  dplyr::arrange(seqnames, start) %>% 
  dplyr::left_join(y = s1Targets, by = unname(tfCols$peakId[bestGrp1Id])) %>% 
  dplyr::left_join(y = s2Targets, by = unname(tfCols$peakId[bestGrp2Id])) %>% 
  dplyr::select(
    seqnames, start, end, name, starts_with("Conc_"), Fold, p.value, FDR, diffBind, peakOccupancy,
    categoryDiffbind, starts_with("peakCall."), starts_with("pvalGood."), contains(bestGrp1Id),
    contains(bestGrp2Id)) %>% 
  dplyr::filter(peakOccupancy != "no_consensus") %>% 
  dplyr::mutate(
    targetMatch = dplyr::case_when(
      !!sym(tfCols$targetGene[bestGrp1Id]) == !!sym(tfCols$targetGene[bestGrp2Id]) ~ TRUE,
      is.na(!!sym(tfCols$targetGene[bestGrp1Id])) & is.na(!!sym(tfCols$targetGene[bestGrp2Id])) ~ TRUE,
      is.na(!!sym(tfCols$targetGene[bestGrp1Id])) ~ NA,
      is.na(!!sym(tfCols$targetGene[bestGrp2Id])) ~ NA,
      categoryDiffbind %in% c(grp1SpecificOcc, grp2SpecificOcc) ~ NA,
      TRUE ~ FALSE
    ),
    peakPosMatch = dplyr::case_when(
      !!sym(tfCols$peakPosition[bestGrp1Id]) == !!sym(tfCols$peakPosition[bestGrp2Id]) ~ TRUE,
      is.na(!!sym(tfCols$peakPosition[bestGrp1Id])) & is.na(!!sym(tfCols$peakPosition[bestGrp2Id])) ~ TRUE,
      is.na(!!sym(tfCols$peakPosition[bestGrp1Id])) ~ NA,
      is.na(!!sym(tfCols$peakPosition[bestGrp2Id])) ~ NA,
      categoryDiffbind %in% c(grp1SpecificOcc, grp2SpecificOcc) ~ NA,
      TRUE ~ FALSE
    )
  )

## get consensus geneId. choose appropriate when there is no consensus target gene
diffAnn <- diffAnn %>% 
  dplyr::mutate(
    geneId = dplyr::case_when(
      peakOccupancy == !!grp1SpecificOcc ~ !! sym(tfCols$targetGene[bestGrp1Id]),
      peakOccupancy == !!grp2SpecificOcc ~ !! sym(tfCols$targetGene[bestGrp2Id]),
      peakOccupancy == "common" & targetMatch ~ !! sym(tfCols$targetGene[bestGrp1Id]), 
      TRUE ~ "NA")
  )


## get consensus peak position. choose appropriate when there is not consensus peakPosition
diffAnn <- diffAnn %>% 
  dplyr::mutate(
    peakPosition = dplyr::case_when(
      peakOccupancy == !!grp1SpecificOcc ~ !! sym(tfCols$peakPosition[bestGrp1Id]),
      peakOccupancy == !!grp2SpecificOcc ~ !! sym(tfCols$peakPosition[bestGrp2Id]),
      peakOccupancy == "common" & peakPosMatch ~ !! sym(tfCols$peakPosition[bestGrp1Id]),
      ## other common peaks
      peakPosMatch & !! sym(tfCols$peakPosition[bestGrp1Id]) == "TES" ~ "TES",
      peakPosMatch & !! sym(tfCols$peakPosition[bestGrp2Id]) == "TES" ~ "TES",
      TRUE ~ "NA")
  )


diffAnn <- diffAnn %>% dplyr::select(seqnames, start, end, name, geneId, peakPosition, everything())

readr::write_tsv(x = diffAnn, file = paste(outPrefix, ".all.annotation.tab", sep = ""))

##################################################################################
## final confident diffbind annotation
## tf1 specific targets
tf1Specific <- dplyr::filter(
  diffAnn,
  !! sym(groupCols$peakCall[grp1]) >= 2 & !! sym(groupCols$peakCall[grp2]) < 2
) %>% 
  dplyr::mutate(
    !! tfCols$hasPeak[bestGrp1Id] := TRUE,
    !! tfCols$hasPeak[bestGrp2Id] := FALSE
  )

## tf2 specific targets
tf2Specific <- dplyr::filter(
  diffAnn,
  !! sym(groupCols$peakCall[grp2]) >= 2 & !! sym(groupCols$peakCall[grp1]) < 2
) %>% 
  dplyr::mutate(
    !! tfCols$hasPeak[bestGrp1Id] := FALSE,
    !! tfCols$hasPeak[bestGrp2Id] := TRUE
  )

## common targets: 
common <- dplyr::filter(
  diffAnn,
  !! sym(groupCols$peakCall[grp2]) >= 2 & !! sym(groupCols$peakCall[grp1]) >= 2,
  targetMatch, peakPosMatch) %>% 
  dplyr::mutate(
    !! tfCols$hasPeak[bestGrp1Id] := TRUE,
    !! tfCols$hasPeak[bestGrp2Id] := TRUE
  )


## no need to select a best peak when a gene has more than 1 peaks. report all
finalDiffbind <- dplyr::bind_rows(tf1Specific, tf2Specific, common) %>% 
  dplyr::arrange(seqnames, start) %>% 
  dplyr::select(-starts_with("targetGene.")) %>% 
  dplyr::select(seqnames, start, end, name, geneId, peakPosition, Fold, p.value, FDR, diffBind,
                peakOccupancy, categoryDiffbind, starts_with("peakCall"), contains(bestGrp1Id),
                contains(bestGrp2Id), everything(), -starts_with("Conc_")
  )


readr::write_tsv(x = finalDiffbind, file = paste(outPrefix, ".annotation.filtered.tab", sep = ""))

##################################################################################





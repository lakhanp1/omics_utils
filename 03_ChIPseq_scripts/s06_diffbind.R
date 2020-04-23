library(chipmine)
library(org.AFumigatus293.eg.db)
library(BSgenome.Afumigatus.AspGD.Af293)
library(DiffBind)

## this script:
## 1) run DiffBind to perform differential binding analysis
## 2) add peak data for individual sample
## 3) prepare combined report
## 4) add peak target gene information for each best peaksets

rm(list = ls())

outDir <- here::here("analysis", "ChIPseq_analysis", "diffBind")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

compare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

analysisName <- "creE_diffbind"
outPrefix <- paste(outDir, "/", analysisName, sep = "")

fdr_cut <- 0.05
lfc_cut <- 1
up_cut <- lfc_cut
down_cut <- lfc_cut * -1
##################################################################################
orgDb <- org.AFumigatus293.eg.db

file_diffbindInfo <- here::here("analysis", "ChIPseq_analysis", "diffBind", "sampleInfo.txt")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")

TF_dataPath <- here::here("data", "TF_data")

##################################################################################
sampleInfo <- suppressMessages(readr::read_tsv(file = file_diffbindInfo))

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleInfo$SampleID,
                                   dataPath = TF_dataPath,
                                   matrixSource = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq", "overlap", "targetGene"),
  FUN = function(x){ structure(paste(x, ".", exptData$sampleId, sep = ""), names = exptData$sampleId) },
  simplify = F, USE.NAMES = T
)

grp1 <- compare[1]
grp2 <- compare[2]
grp1Index <- which(exptData$groupId == grp1)
grp2Index <- which(exptData$groupId == grp2)
grp1Samples <- exptData$sampleId[grp1Index]
grp2Samples <- exptData$sampleId[grp2Index]
grp1SpecificOcc = paste("specific:", grp1, sep = "")
grp2SpecificOcc = paste("specific:", grp2, sep = "")

groupCols <- sapply(
  X = c("peakCall"),
  FUN = function(x){ structure(paste(x, ".", compare, sep = ""), names = compare) },
  simplify = F, USE.NAMES = T
)

##################################################################################
## diffBind analysis

# dataDba <- DiffBind::dba(sampleSheet = sampleInfo)
# 
# plot(dataDba)
# 
# countsDba <- DiffBind::dba.count(DBA = dataDba, score = "DBA_SCORE_TMM_MINUS_FULL")
# 
# plot(countsDba)
# 
# countsDba <- DiffBind::dba.contrast(DBA = countsDba, categories = DBA_CONDITION, minMembers = 2)
# 
# diffDba <- DiffBind::dba.analyze(DBA = countsDba, method = DBA_DESEQ2, bReduceObjects = FALSE)
# 
# DiffBind::dba.save(DBA = diffDba, file = paste(analysisName, ".dba", sep = ""), dir = outDir, pre = "")

## load DBA object
diffDba <- DiffBind::dba.load(file = paste(analysisName, ".dba", sep = ""), dir = outDir, pre = "")

plot(diffDba, contrast=1)

diffDf <- DiffBind::dba.report(DBA = diffDba, bFlip = TRUE, th = 1,
                               bCalledDetail = T, DataType = DBA_DATA_FRAME)

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
diffGr <- DiffBind::dba.report(DBA = diffDba, bFlip = TRUE, th = 1)
mcols(diffGr)$name <- paste("peak_region", 1:length(diffGr), sep = "_")

rtracklayer::export.bed(object = diffGr, con = paste(outPrefix, ".merged_regions.bed", sep = ""),
                        format = "bed")

## add overlapping peaks from each sample
diffRes <- combinatorial_binding_matrix(sampleInfo = exptData, peakRegions = diffGr)

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
      !!as.name(groupCols$peakCall[grp1]) >= 2 & !!as.name(groupCols$peakCall[grp2]) >= 2 ~ "common",
      !!as.name(groupCols$peakCall[grp1]) >= 2 & !!as.name(groupCols$peakCall[grp2]) < 2 ~ 
        !! grp1SpecificOcc,
      !!as.name(groupCols$peakCall[grp1]) < 2 & !!as.name(groupCols$peakCall[grp2]) >= 2 ~ 
        !! grp2SpecificOcc,
      TRUE ~ "no_consensus"
    )
  )


dplyr::group_by(diffData, peakOccupancy, diffBind) %>% 
  dplyr::summarise(n = n())

## assign peak category based on diffbind and peakOccupancy
diffData <- diffData %>% 
  dplyr::mutate(
    categoryDiffbind = dplyr::case_when(
      diffBind == "down" | peakOccupancy == grp1SpecificOcc ~ paste(grp1, ":enriched", sep = ""),
      diffBind == "up" | peakOccupancy == grp2SpecificOcc ~ paste(grp2, ":enriched", sep = ""),
      diffBind == "noDiff" & peakOccupancy == "common" ~ "common",
      TRUE ~ "NA"
    )
  ) %>% 
  dplyr::select(seqnames, start, end, name, starts_with("Conc"), Fold, p.value, FDR,
                diffBind, peakOccupancy, categoryDiffbind, starts_with("peakCall."), everything(),
                -starts_with("peakEnrichment."))


readr::write_tsv(x = diffData, path = paste(outPrefix, ".all.tab", sep = ""))

##################################################################################
## diffbind report with target genes

bestGrp1Id <- exptData$sampleId[exptData$groupId == grp1 & exptData$bestRep == 1]
bestGrp2Id <- exptData$sampleId[exptData$groupId == grp2 & exptData$bestRep == 1]

## if a condition has no consensus peaks in both the replicates, 
## remove the peak information even if best replicate has peak
diffData <- diffData %>% 
  dplyr::mutate(
    !!as.name(tfCols$peakId[bestGrp1Id]) := ifelse(
      !!as.name(groupCols$peakCall[grp1]) == 1, yes = NA, no = !!as.name(tfCols$peakId[bestGrp1Id])
    ),
    !!as.name(tfCols$overlap[bestGrp1Id]) := ifelse(
      !!as.name(groupCols$peakCall[grp1]) == 1, yes = FALSE, no = !!as.name(tfCols$overlap[bestGrp1Id])
    ),
    !!as.name(tfCols$peakPval[bestGrp1Id]) := ifelse(
      !!as.name(groupCols$peakCall[grp1]) == 1, yes = NA, no = !!as.name(tfCols$peakPval[bestGrp1Id])
    )
  ) %>%
  dplyr::mutate(
    !!as.name(tfCols$peakId[bestGrp2Id]) := ifelse(
      !!as.name(groupCols$peakCall[grp2]) == 1, yes = NA, no = !!as.name(tfCols$peakId[bestGrp2Id])
    ),
    !!as.name(tfCols$overlap[bestGrp2Id]) := ifelse(
      !!as.name(groupCols$peakCall[grp2]) == 1, yes = FALSE, no = !!as.name(tfCols$overlap[bestGrp2Id])
    ),
    !!as.name(tfCols$peakPval[bestGrp2Id]) := ifelse(
      !!as.name(groupCols$peakCall[grp2]) == 1, yes = NA, no = !!as.name(tfCols$peakPval[bestGrp2Id])
    )
  ) 

## add pvalFiltered column
diffData <- diffData %>% 
  dplyr::mutate(
    !! as.name(tfCols$pvalFiltered[bestGrp1Id]) := if_else(
      condition = !! as.name(tfCols$peakPval[bestGrp1Id]) >= exptDataList[[bestGrp1Id]]$pval_cutoff,
      true = TRUE, false = FALSE
    ),
    !! as.name(tfCols$pvalFiltered[bestGrp2Id]) := if_else(
      condition = !! as.name(tfCols$peakPval[bestGrp2Id]) >= exptDataList[[bestGrp2Id]]$pval_cutoff,
      true = TRUE, false = FALSE
    )
  )


## import peak annotation for sample1 and add new column with pval cutoff pass result
s1Targets <- import_peak_annotation(
  sampleId = bestGrp1Id,
  peakAnnoFile = exptDataList[[bestGrp1Id]]$narrowpeakAnno,
  columns = c("peakId", "geneId", "peakPosition", "peakType", "peakDist")
) %>% 
  dplyr::filter(!is.na(geneId)) %>% 
  dplyr::filter(!grepl(pattern = "pseudo_", x = !!as.name(tfCols$peakType[bestGrp1Id]))) %>% 
  dplyr::rename(!! as.name(tfCols$targetGene[bestGrp1Id]) := geneId)


## import peak annotation for sample2 and add new column with pval cutoff pass result
s2Targets <- import_peak_annotation(
  sampleId = bestGrp2Id,
  peakAnnoFile = exptDataList[[bestGrp2Id]]$narrowpeakAnno,
  columns = c("peakId", "geneId", "peakPosition", "peakType", "peakDist")
) %>% 
  dplyr::filter(!is.na(geneId)) %>% 
  dplyr::filter(!grepl(pattern = "pseudo_", x = !!as.name(tfCols$peakType[bestGrp2Id]))) %>% 
  dplyr::rename(!! as.name(tfCols$targetGene[bestGrp2Id]) := geneId)



## combine diffbind results with peak target annotation information
diffAnn <- diffData %>% 
  dplyr::arrange(seqnames, start) %>% 
  dplyr::left_join(y = s1Targets, by = unname(tfCols$peakId[bestGrp1Id])) %>% 
  dplyr::left_join(y = s2Targets, by = unname(tfCols$peakId[bestGrp2Id])) %>% 
  dplyr::select(seqnames, start, end, name, starts_with("Conc_"), Fold, p.value, FDR, diffBind, peakOccupancy,
                categoryDiffbind, starts_with("peakCall."), contains(bestGrp1Id), contains(bestGrp2Id)) %>% 
  dplyr::filter(peakOccupancy != "no_consensus") %>% 
  dplyr::mutate(
    targetMatch = if_else(
      condition = !! as.name(tfCols$targetGene[bestGrp1Id]) == !! as.name(tfCols$targetGene[bestGrp2Id]),
      true = TRUE, false = FALSE),
    peakPosMatch = if_else(
      condition = !! as.name(tfCols$peakPosition[bestGrp1Id]) == !! as.name(tfCols$peakPosition[bestGrp2Id]),
      true = TRUE, false = FALSE)
  )

## get consensus geneId. choose appropriate when there is no consensus target gene
diffAnn <- diffAnn %>% 
  dplyr::mutate(
    geneId = dplyr::case_when(
      peakOccupancy == !!grp1SpecificOcc ~ !! as.name(tfCols$targetGene[bestGrp1Id]),
      peakOccupancy == !!grp2SpecificOcc ~ !! as.name(tfCols$targetGene[bestGrp2Id]),
      peakOccupancy == "common" & targetMatch ~ !! as.name(tfCols$targetGene[bestGrp1Id]), 
      TRUE ~ "NA")
  )

## get consensus peak position. choose appropriate when there is not consensus peakPosition
diffAnn <- diffAnn %>% 
  dplyr::mutate(
    peakPosition = dplyr::case_when(
      peakOccupancy == !!grp1SpecificOcc ~ !! as.name(tfCols$peakPosition[bestGrp1Id]),
      peakOccupancy == !!grp2SpecificOcc ~ !! as.name(tfCols$peakPosition[bestGrp2Id]),
      peakOccupancy == "common" & peakPosMatch ~ !! as.name(tfCols$peakPosition[bestGrp1Id]),
      ## other common peaks
      peakPosMatch & !! as.name(tfCols$peakPosition[bestGrp1Id]) == "TES" ~ "TES",
      peakPosMatch & !! as.name(tfCols$peakPosition[bestGrp2Id]) == "TES" ~ "TES",
      TRUE ~ "NA")
  )


## get bestPval between two samples. this will be used to select best peak pair when a gene has
## multiple peaks
diffAnn <- diffAnn %>% 
  dplyr::mutate(
    bestPval = dplyr::case_when(
      peakOccupancy == !!grp1SpecificOcc ~ !!as.name(tfCols$peakPval[bestGrp1Id]),
      peakOccupancy == !!grp2SpecificOcc ~ !!as.name(tfCols$peakPval[bestGrp2Id]),
      peakOccupancy == "common" &
        !!as.name(tfCols$peakPval[bestGrp1Id]) > !!as.name(tfCols$peakPval[bestGrp2Id]) ~
        !!as.name(tfCols$peakPval[bestGrp1Id]),
      peakOccupancy == "common" &
        !!as.name(tfCols$peakPval[bestGrp1Id]) < !!as.name(tfCols$peakPval[bestGrp2Id]) ~
        !!as.name(tfCols$peakPval[bestGrp2Id]),
      TRUE ~ 0)
  )


diffAnn <- diffAnn %>% dplyr::select(seqnames, start, end, name, geneId, peakPosition, everything())

diffAnn$pvalFilteredN <- purrr::pmap_int(
  .l = dplyr::select(diffAnn, starts_with("pvalFiltered.")),
  .f = sum, na.rm = TRUE
)

readr::write_tsv(x = diffAnn, path = paste(outPrefix, ".all.annotation.tab", sep = ""))

##################################################################################
## final confident diffbind annotation
## tf1 specific targets. additionally remove the tf1 specific targets which are
## no-diff by diffbind and there is one peak detected in another sample
tf1Specific <- dplyr::filter(
  diffAnn,
  !! as.name(groupCols$peakCall[grp1]) == 2 & !! as.name(groupCols$peakCall[grp2]) < 2,
  !is.na(!! as.name(tfCols$peakPosition[bestGrp1Id]))) %>% 
  dplyr::filter(
    !(diffBind == "noDiff" & !! as.name(groupCols$peakCall[grp2]) == 1)
  ) %>% 
  dplyr::mutate(
    !! tfCols$hasPeak[bestGrp1Id] := TRUE,
    !! tfCols$hasPeak[bestGrp2Id] := FALSE
  )

## tf2 specific targets. additionally remove the tf2 specific targets which are
## no-diff by diffbind and there is one peak detected in another sample
tf2Specific <- dplyr::filter(
  diffAnn,
  !! as.name(groupCols$peakCall[grp2]) == 2 & !! as.name(groupCols$peakCall[grp1]) < 2,
  !is.na(!! as.name(tfCols$peakPosition[bestGrp2Id]))) %>% 
  dplyr::filter(
    !(diffBind == "noDiff" & !! as.name(groupCols$peakCall[grp1]) == 1)
  ) %>% 
  dplyr::mutate(
    !! tfCols$hasPeak[bestGrp1Id] := FALSE,
    !! tfCols$hasPeak[bestGrp2Id] := TRUE
  )

## common targets
common <- dplyr::filter(
  diffAnn,
  !! as.name(groupCols$peakCall[grp2]) == 2 & !! as.name(groupCols$peakCall[grp1]) == 2,
  targetMatch, peakPosMatch) %>% 
  dplyr::mutate(
    !! tfCols$hasPeak[bestGrp1Id] := TRUE,
    !! tfCols$hasPeak[bestGrp2Id] := TRUE
  )


finalDiffbind <- dplyr::bind_rows(tf1Specific, tf2Specific, common) %>% 
  dplyr::group_by(geneId, peakPosition) %>% 
  dplyr::arrange(desc(bestPval),  .by_group = TRUE) %>%
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(seqnames, start) %>% 
  dplyr::select(-starts_with("targetGene.")) %>% 
  dplyr::select(seqnames, start, end, name, geneId, peakPosition, Fold, p.value, FDR, diffBind,
                peakOccupancy, categoryDiffbind, starts_with("peakCall"), contains(bestGrp1Id),
                contains(bestGrp2Id), everything(), -starts_with("Conc_")
  )


readr::write_tsv(x = finalDiffbind, path = paste(outPrefix, ".annotation.filtered.tab", sep = ""))

##################################################################################





library(chipmine)
library(here)
library(BSgenome.Afumigatus.AspGD.Af293)
library(data.table)
library(org.AFumigatus293.eg.db)
library(TxDb.Afumigatus.Af293.AspGD.GFF)

## this script
## generates peak target matrix from diffbind annotation results
rm(list = ls())

outDir <- here::here("analysis", "02_ChIPseq_analysis", "01_peak_targets")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, "/", sep = "")

compare <- c("CREEHA_CONTROL", "CREEHA_10MMAA")

##################################################################################
summitSeqLen <- 500

file_diffbindInfo <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind", "sampleInfo.txt")
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_diffbindRes <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind", "creE_diffbind.all.annotation.tab")
file_diffbindAnn <- here::here("analysis", "02_ChIPseq_analysis", "03_diffBind", "creE_diffbind.annotation.filtered.tab")
TF_dataPath <- here::here("data", "TF_data")

sampleInfo <- suppressMessages(readr::read_tsv(file = file_diffbindInfo))

orgDb <- org.AFumigatus293.eg.db
txDb <- TxDb.Afumigatus.Af293.AspGD.GFF

## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   dataPath = TF_dataPath,
                                   profileMatrixSuffix = "normalizedmatrix")

exptDataList <- purrr::transpose(exptData) %>%
  purrr::set_names(nm = purrr::map(., "sampleId"))

grp1 <- compare[1]
grp2 <- compare[2]

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq", "summitRegion"),
  FUN = function(x){ structure(paste(x, ".", exptData$sampleId, sep = ""), names = exptData$sampleId) },
  simplify = F, USE.NAMES = T
)

geneStrands <- as.data.frame(genes(txDb), row.names = NULL) %>% 
  dplyr::select(geneId = gene_id, strand)

geneInfo <- AnnotationDbi::select(x = orgDb, keys = keys(orgDb), columns = "DESCRIPTION") %>% 
  dplyr::rename(geneId = GID) %>% 
  dplyr::left_join(geneStrands, by = "geneId") %>% 
  dplyr::select(geneId, strand, DESCRIPTION)



##################################################################################

diffbindAnn <- suppressMessages(readr::read_tsv(file = file_diffbindAnn, col_names = T)) %>% 
  dplyr::mutate(DiffBind_region = paste(seqnames, ":", start, "-", end, sep = "")) %>% 
  dplyr::rename(DiffBind.foldChange = Fold,
                DiffBind.pvalue = p.value,
                DiffBind.FDR = FDR)

# diffbindAnn <- dplyr::filter(diffbindAnn, pvalFilteredN > 0)

## prepare ChIPseq target data
bestGrp1Id <- exptData$sampleId[exptData$groupId == grp1 & exptData$bestRep == 1]
bestGrp2Id <- exptData$sampleId[exptData$groupId == grp2 & exptData$bestRep == 1]

tf1SummitSeq <- get_narrowpeak_summit_seq(npFile = exptDataList[[bestGrp1Id]]$peakFile,
                                          id = exptDataList[[bestGrp1Id]]$sampleId,
                                          genome = BSgenome.Afumigatus.AspGD.Af293,
                                          length = summitSeqLen)

tf2SummitSeq <- get_narrowpeak_summit_seq(npFile = exptDataList[[bestGrp2Id]]$peakFile,
                                          id = exptDataList[[bestGrp2Id]]$sampleId,
                                          genome = BSgenome.Afumigatus.AspGD.Af293,
                                          length = summitSeqLen)

targetMat <- diffbindAnn %>% 
  dplyr::left_join(y = tf1SummitSeq, by = structure("name", names = tfCols$peakId[bestGrp1Id])) %>%
  dplyr::left_join(y = tf2SummitSeq, by = structure("name", names = tfCols$peakId[bestGrp2Id])) %>%
  dplyr::select(geneId, name, DiffBind_region, DiffBind.foldChange, DiffBind.pvalue, DiffBind.FDR,
                peakPosition, diffBind, peakOccupancy, categoryDiffbind, pvalFilteredN,
                starts_with("peakCall."),
                starts_with("hasPeak."), starts_with("peakPval."), starts_with("peakId."),
                starts_with("summitSeq."), starts_with("summitRegion."),
                starts_with("peakType."), starts_with("peakDist."), bestPval) %>% 
  dplyr::select(geneId, name, DiffBind_region, DiffBind.foldChange, DiffBind.pvalue, DiffBind.FDR,
                peakPosition, diffBind, peakOccupancy, categoryDiffbind, starts_with("peakCall."),
                contains(bestGrp1Id), contains(bestGrp2Id), pvalFilteredN, bestPval)

targetMat <- dplyr::left_join(x = targetMat, y = geneInfo, by = "geneId")

readr::write_tsv(x = targetMat, path = paste(outPrefix, "diffbind_allPeak_targets.tab", sep = ""))

##################################################################################
## get individual peak target gene
peakAnnotation <- suppressMessages(readr::read_tsv(file = file_diffbindRes, col_names = T)) %>% 
  dplyr::select(starts_with("peakId."), geneId, peakPosition) %>% 
  dplyr::filter(!is.na(geneId)) %>% 
  as.data.table() %>% 
  data.table::melt.data.table(
    measure.vars = purrr::map(.x = tfCols, .f = function(x) x[c(bestGrp1Id, bestGrp2Id)]) %>% 
      .[c("peakId")],
    variable.name = "sampleId", value.name = "peakId"
  ) %>% 
  tibble::as_tibble() %>% 
  dplyr::filter(!is.na(peakId)) %>% 
  dplyr::select(-sampleId)

## get diffbind results
diffbindRes <- suppressMessages(readr::read_tsv(file = file_diffbindRes, col_names = T)) %>% 
  dplyr::select(name, diffBind, peakOccupancy, categoryDiffbind, pvalFilteredN,
                starts_with("peakId."), starts_with("peakPval."), starts_with("pvalFiltered.")) %>% 
  dplyr::left_join(y = tf1SummitSeq, by = structure("name", names = tfCols$peakId[bestGrp1Id])) %>%
  dplyr::left_join(y = tf2SummitSeq, by = structure("name", names = tfCols$peakId[bestGrp2Id])) %>% 
  dplyr::distinct()

longTargetDf <- data.table::melt.data.table(
  data = as.data.table(diffbindRes),
  measure.vars = purrr::map(.x = tfCols, .f = function(x) x[c(bestGrp1Id, bestGrp2Id)]) %>% 
    .[c("peakId", "peakPval", "pvalFiltered", "summitRegion", "summitSeq")],
  variable.name = "sampleId"
)

levels(longTargetDf$sampleId) <- c(bestGrp1Id, bestGrp2Id)

peakwiseDiff <- as_tibble(longTargetDf) %>% 
  dplyr::filter(!is.na(peakId)) %>% 
  dplyr::left_join(y = peakAnnotation, by = "peakId") %>% 
  dplyr::left_join(y = geneStrands, by = "geneId") %>% 
  dplyr::select(peakId, peakPval, pvalFiltered, sampleId, everything()) %>% 
  dplyr::distinct()

summitSeq <- DNAStringSet(x = peakwiseDiff$summitSeq)
negStrandTargets <- which(peakwiseDiff$strand == "-")
summitSeq[negStrandTargets] <- reverseComplement(summitSeq[negStrandTargets])
peakwiseDiff$summitSeqStranded <- as.character(summitSeq, use.names = FALSE)

readr::write_tsv(x = peakwiseDiff, path = paste(outPrefix, "peakwise_diffbind_groups.tab", sep = ""))


motifData <- dplyr::group_by(peakwiseDiff, name) %>% 
  dplyr::arrange(desc(peakPval)) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup()

readr::write_tsv(x = motifData, path = paste(outPrefix, "summit_seq_data.tab", sep = ""))



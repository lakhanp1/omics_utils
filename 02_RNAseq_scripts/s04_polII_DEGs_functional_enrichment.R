suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(require(openxlsx))
suppressPackageStartupMessages(library(argparse))


## polII ChIPseq DEG functional analysis using topGO, clusterProfiler and fgsea

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

diffDataPath <- here::here("analysis", "06_polII_diff")

cutoff_fpkm <- 10
cutoff_fdr <- 0.05
cutoff_lfc <- 1
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

col_lfc <- "log2FoldChange"

orgDb <- org.Anidulans.FGSCA4.eg.db
keggOrg <- 'ani'
col_degId <- "geneId"
col_degOrgdbKey <- "GID"
col_kegg <- "KEGG_ID"
col_gsea <- "geneId"
col_geneName <- "GENE_NAME"
# file_msigDesc <- "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/msigDB_geneset_desc.tab"

###########################################################################
#############################################
## Run DESeq2 pipeline using Rscript       ##
## command line arguments                  ##
#############################################

parser <- ArgumentParser(
  description = "This script automates the functional enrichment analysis of RNAseq DEG results."
)

parser$add_argument(
  "-c", "--config", action="store",
  dest = "config", type = "character", nargs = 1, required = TRUE,
  # default = here::here("data", "reference_data", "polII_DESeq2_info.txt"),
  help = "DEG configuration TAB delimited file with columns: comparison, type, group1, group2, design, samples"
)

parser$add_argument(
  "-d", "--deg", action="store",
  dest = "deg", type = "character", nargs = 1, required = TRUE,
  help = "DEG comparison ID. This value should be present in the column comparison of config file"
)

# parser$print_help()

args <- parser$parse_args()

file_RNAseq_info <- args$config
degResult <- args$deg

# file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
# degResult <- "AN10295_sCopy_OE_vs_WT"

outDir <- paste(diffDataPath, "/", degResult, sep = "")
outPrefix <- paste(outDir, "/", degResult, sep = "")

###########################################################################

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == degResult)

if(nrow(rnaseqInfo) != 1){
  stop(analysisName, " RNAseq data does not exist in config file")
}

file_degs <- paste(outPrefix, ".DEG_FPKM_filtered.txt", sep = "")

degs <- suppressMessages(readr::read_tsv(file = file_degs)) %>% 
  dplyr::mutate(
    fpkmFilterSign = dplyr::if_else(
      condition = maxFpkm >= cutoff_fpkm, true = 1 * sign(shrinkLog2FC),
      false = 0, missing = 0
    ),
    rankMetric = (-log10(pvalue) * sign(shrinkLog2FC))
  ) %>%
  dplyr::arrange(desc(rankMetric)) %>% 
  dplyr::filter(!is.na(rankMetric))

missingCols <- setdiff(c(col_gsea, col_kegg), colnames(degs))

if(!is.null(missingCols)){
  gseaInfo <- suppressMessages(
    AnnotationDbi::select(
      x = orgDb, keys = degs$geneId, keytype = col_degOrgdbKey,
      columns = missingCols
    )
  ) %>% 
    dplyr::filter_at(.vars = missingCols, .vars_predicate = all_vars(!is.na(.))) %>% 
    dplyr::rename(geneId = !!col_degOrgdbKey)
  
  degs <- dplyr::left_join(x = degs, y = gseaInfo, by = "geneId")
}

downDegs <- dplyr::filter(
  degs, padj <= cutoff_fdr & !!sym(col_lfc) <= cutoff_down, maxFpkm >= cutoff_fpkm
) %>% 
  dplyr::mutate(category = "down")

upDegs <- dplyr::filter(
  degs, padj <= cutoff_fdr & !!sym(col_lfc) >= cutoff_up, maxFpkm >= cutoff_fpkm
) %>% 
  dplyr::mutate(category = "up")

degData <- dplyr::bind_rows(upDegs, downDegs)

contrast <- unique(upDegs$comparison)

geneList <- dplyr::filter(degs, !is.na(!!sym(col_gsea))) %>% 
  dplyr::select(!!col_gsea, rankMetric) %>% 
  tibble::deframe()

## replace +Inf and -Inf values with max and min
geneList[is.infinite(geneList) & geneList > 0] <- max(geneList[is.finite(geneList)]) + 1

geneList[is.infinite(geneList) & geneList < 0] <- min(geneList[is.finite(geneList)]) - 1

###########################################################################
## topGO GO enrichment
topgo_up <- topGO_enrichment(
  genes = unique(upDegs$geneId),
  orgdb = orgDb, inKeytype = col_degOrgdbKey,
  type = "BP", goNodeSize = 5,
  genenameKeytype = col_geneName
)

topgo_down <- topGO_enrichment(
  genes = unique(downDegs$geneId),
  orgdb = orgDb, inKeytype = col_degOrgdbKey,
  type = "BP", goNodeSize = 5,
  genenameKeytype = col_geneName
)

topgo_res <- dplyr::bind_rows(
  if(!is.null(topgo_up)) dplyr::mutate(.data = as.data.frame(topgo_up), category = "up"),
  if(!is.null(topgo_down)) dplyr::mutate(.data = as.data.frame(topgo_down), category = "down")
)

# dplyr::glimpse(topgo_res)
out_topGO <- paste(outPrefix, ".topGO.tab", sep = "")

if(nrow(topgo_res) != 0){
  topgo_res <- dplyr::mutate(topgo_res, contrast = contrast)
  readr::write_tsv(x = topgo_res, file = out_topGO)
} else{
  readr::write_file(x = topgo_res, file = out_topGO)
}

# ## top 10 GO term bar plot
# topgoPlotDf <- dplyr::group_by(topgo_res, category) %>%
#   dplyr::arrange(pvalue, .by_group = TRUE) %>%
#   dplyr::slice(1:10) %>%
#   dplyr::ungroup()
#
#
# topgo_bar <- enrichment_bar(
#   df = topgoPlotDf,
#   title = paste(degResult, "\ntop 10 enriched GO terms in up and down DEGs")
# )
#
# png(filename = paste(outPrefix, ".topGO_bar.png", sep = ""),
#     width = 2500, height = 2500, res = 250)
#
# topgo_bar
# dev.off()

###########################################################################

kegg_up <- fgsea_kegg_overrepresentation(
  genes = unique(upDegs$geneId), keggOrg = keggOrg, orgdb = orgDb,
  inKeytype = col_degOrgdbKey, keggKeytype = col_kegg,
  genenameKeytype = col_geneName, pvalueCutoff = 0.05
)

kegg_down <- fgsea_kegg_overrepresentation(
  genes = unique(downDegs$geneId), keggOrg = keggOrg, orgdb = orgDb,
  inKeytype = col_degOrgdbKey, keggKeytype = col_kegg,
  genenameKeytype = col_geneName, pvalueCutoff = 0.05
)


kegg_res <- dplyr::bind_rows(
  if(!is.null(kegg_up)) dplyr::mutate(.data = as.data.frame(kegg_up), category = "up"),
  if(!is.null(kegg_down)) dplyr::mutate(.data = as.data.frame(kegg_down), category = "down")
)

out_kegg <- paste(outPrefix, ".KEGG_fora.tab", sep = "")

if(nrow(kegg_res) != 0){
  kegg_res <- dplyr::mutate(kegg_res, contrast = contrast)

  ## use data.table::fwrite because of list columns
  data.table::fwrite(
    x = kegg_res, file = out_kegg,
    sep = "\t", sep2 = c("",";",""), eol = "\n"
  )
} else{
  readr::write_file(x = kegg_res, file = out_kegg)
}


# ## top 10 KEGG pathway bar plot
# keggPlotDf <- dplyr::group_by(kegg_res, category) %>%
#   dplyr::arrange(pvalue, .by_group = TRUE) %>%
#   dplyr::slice(1:10) %>%
#   dplyr::ungroup()
#
#
# kegg_bar <- enrichment_bar(
#   df = keggPlotDf,
#   title = paste(degResult, "\ntop 10 enriched KEGG pathways in up and down DEGs"),
#   pvalCol = "pval", termCol = "description",
#   colorCol = "category", countCol = "Significant"
# )
#
# png(filename = paste(outPrefix, ".KEGG_bar.png", sep = ""),
#     width = 2500, height = 2500, res = 250)
#
# kegg_bar
# dev.off()

###########################################################################
## GSEA enrichment using fgsea
goGsea <- fgsea_orgdb_GO_KEGG(
  genelist = geneList, orgdb = orgDb, inKeytype = col_degOrgdbKey,
  keggOrg = keggOrg, keggKeytype = col_kegg, genenameKeytype = col_geneName,
  minNodeSize = 10, maxNodeSize = 1000, eps = 0
)

goGsea_res <- dplyr::mutate(goGsea, contrast = contrast) %>%
  dplyr::arrange(pval) %>%
  dplyr::select(pathway, description, pathway_type, everything(), -starts_with("leadingEdge"),
                starts_with("leadingEdge"), -contrast, contrast)

# topPathways <- uniqueFgsea[head(order(pval), n=20)][order(NES), pathway]
#
# pt_gsea <- plotGseaTable(
#   goGsea$goGeneset[topPathways], geneList,
#   uniqueFgsea, gseaParam=0.5,
#   colwidths = c(10, 5, 1, 1.5, 1.5),
#   render = FALSE
# )
#
# png(filename = paste(outPrefix, ".fgsea_plot.png", sep = ""),
#     width = 2500, height = 2500, res = 250)
# grid::grid.draw(pt_gsea)
# dev.off()

out_gsea <- paste(outPrefix, ".GO_KEGG_fgsea.tab", sep = "")

## use data.table::fwrite because of list columns
if(nrow(goGsea_res) != 0){
  data.table::fwrite(
    x = goGsea_res, file = out_gsea,
    sep = "\t", sep2 = c("",";",""), eol = "\n"
  )
} else{
  readr::write_file(x = kegg_res, file = out_gsea)
}


###########################################################################
## write data to excel file
descString <- paste(
  "## log2FoldChange cutoff =", cutoff_up, "(up) /", cutoff_down, "(down)",
  "; FPKM cutoff =", cutoff_fpkm, "; q-value cutoff =", cutoff_fdr)

wb <- openxlsx::createWorkbook(creator = "Lakhansing Pardeshi Genomics Core")
headerStyle <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#e6e6e6")

# topgo_res <- suppressMessages(readr::read_tsv(file = paste(outPrefix, ".topGO.tab", sep = "")))
# kegg_res <- suppressMessages(readr::read_tsv(file = paste(outPrefix, ".KEGG_fora.tab", sep = "")))
# goGsea_res <- suppressMessages(readr::read_tsv(file = paste(outPrefix, ".GO_KEGG_fgsea.tab", sep = "")))

wrkSheet <- "topGO"
openxlsx::addWorksheet(wb = wb, sheetName = wrkSheet)
openxlsx::writeData(
  wb = wb, sheet = wrkSheet, startCol = 2, startRow = 1,
  x = paste("## GO enrichment using topGO:", degResult, descString)
)
openxlsx::writeData(
  wb = wb, sheet = wrkSheet, x = topgo_res,
  startCol = 1, startRow = 2, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::addStyle(wb = wb, sheet = wrkSheet, style = headerStyle, rows = 2, cols = 1:ncol(topgo_res))
openxlsx::setColWidths(wb = wb, sheet = wrkSheet, cols = 1, widths = "auto")
openxlsx::setColWidths(wb = wb, sheet = wrkSheet, cols = 2, widths = 60)
openxlsx::freezePane(wb = wb, sheet = wrkSheet, firstActiveRow = 3, firstActiveCol = 2)


wrkSheet <- "KEGG"
openxlsx::addWorksheet(wb = wb, sheetName = wrkSheet)
openxlsx::writeData(
  wb = wb, sheet = wrkSheet, startCol = 2, startRow = 1,
  x = paste("KEGG pathway enrichment using fgsea::fora:", degResult, descString)
)
openxlsx::writeData(
  wb = wb, sheet = wrkSheet, x = kegg_res,
  startCol = 1, startRow = 2, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::addStyle(wb = wb, sheet = wrkSheet, style = headerStyle, rows = 2, cols = 1:ncol(kegg_res))
openxlsx::setColWidths(wb = wb, sheet = wrkSheet, cols = 1, widths = "auto")
openxlsx::setColWidths(wb = wb, sheet = wrkSheet, cols = 2, widths = 60)
openxlsx::freezePane(wb = wb, sheet = wrkSheet, firstActiveRow = 3, firstActiveCol = 2)


wrkSheet <- "fgsea"
openxlsx::addWorksheet(wb = wb, sheetName = wrkSheet)
openxlsx::writeData(
  wb = wb, sheet = wrkSheet, startCol = 2, startRow = 1,
  x = paste("preranked gene set enrichment analysis (GSEA) of GO terms and KEGG pathways:", degResult)
)
openxlsx::writeData(
  wb = wb, sheet = wrkSheet, x = goGsea_res,
  startCol = 1, startRow = 2, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::addStyle(wb = wb, sheet = wrkSheet, style = headerStyle, rows = 2, cols = 1:ncol(goGsea_res))
openxlsx::setColWidths(wb = wb, sheet = wrkSheet, cols = 1, widths = "auto")
openxlsx::setColWidths(wb = wb, sheet = wrkSheet, cols = 2, widths = 60)
openxlsx::freezePane(wb = wb, sheet = wrkSheet, firstActiveRow = 3, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(wb = wb, file = paste(outPrefix, ".enrichment.xlsx", sep = ""), overwrite = TRUE)

###########################################################################




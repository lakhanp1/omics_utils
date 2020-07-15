suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.HSapiens.gencodev30.eg.db))
suppressPackageStartupMessages(library(argparse))


##
## This script plots the 
## 1) compare the two DEG list and plot Upset plot
## 2) perform GO enrichment for each subset of the genes

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")
source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/s01_topGO_functions.R")

###########################################################################

file_RNAseq_info <- here::here("data", "reference_data", "DESeq2_DEG_info.txt")
diffDataPath <- here::here("analysis", "02_DESeq2_diff")
degCompPath <- here::here("analysis", "04_DEG_compare")

cutoff_fdr <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

col_lfc <- "log2FoldChange"

orgDb <- org.HSapiens.gencodev30.eg.db
keggOrg <- 'hsa'
col_degOrgdbKey <- "ENSEMBL_VERSION"
col_kegg <- "NCBI_ID"
col_gsea <- "NCBI_ID"
col_topGO <- "ENSEMBL"
col_geneName <- "GENE_NAME"
file_topGO <- "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/geneid2go.HSapiens.GRCh38p12.topGO.map"
file_msigDesc <- "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/msigDB_geneset_desc.tab"

###########################################################################
#############################################
## Run DESeq2 pipeline using Rscript       ##
## command line arguments                  ##
#############################################

parser <- ArgumentParser(
  description = "This script automates the comparison analysis of two DEG list."
)

parser$add_argument(
  "-c", "--config", action="store",
  dest = "config", type = "character", nargs = 1, required = TRUE,
  # default = here::here("data", "reference_data", "polII_DESeq2_info.txt"),
  help = "TAB delimited configuration file with columns: degPairId, deg1, deg2"
)

parser$add_argument(
  "-p", "--pair", action="store",
  dest = "pair", type = "character", nargs = 1, required = TRUE,
  help = "DEG pair comparison ID. This value should be present in the degPairId column of config file"
)

# parser$print_help()

args <- parser$parse_args()

file_config <- args$config
analysisName <- args$pair

# file_config <- here::here("analysis", "04_DEG_compare", "DEG_pair_compare.conf.tab")
# analysisName <- "5uM_12hr"

cat("DEG pair ID: ", analysisName, "\n")

###########################################################################
paircompConf <- suppressMessages(readr::read_tsv(file = file_config)) %>%
  dplyr::filter(degPairId == analysisName)

degResults <-  c(paircompConf$deg1, paircompConf$deg2)

outDir <- paste(degCompPath, analysisName, sep = "/")
outPrefix <- paste(outDir, analysisName, sep = "/")

if(!dir.exists(outDir)){
  dir.create(outDir)
}

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degResults)


## use org.db
geneInfo <- AnnotationDbi::select(x = orgDb,
                                  keys = keys(x = orgDb, keytype = "ENSEMBL_VERSION"),
                                  columns = c("ENSEMBL", "GENE_NAME", "DESCRIPTION"),
                                  keytype = "ENSEMBL_VERSION") %>% 
  dplyr::rename(geneId = ENSEMBL_VERSION)


###########################################################################
## import data

## function to extract the log2FoldChange, padj and diff coulumns for each DEG result file
get_foldchange <- function(degFile, name, lfcCol = "log2FoldChange", otherCols = NULL){
  
  degs <- suppressMessages(readr::read_tsv(file = degFile)) %>% 
    dplyr::mutate(
      diff = dplyr::case_when(
        !!sym(lfcCol) >= cutoff_up & padj <= cutoff_fdr ~ "up",
        !!sym(lfcCol) <= cutoff_down & padj <= cutoff_fdr ~ "down",
        TRUE ~ "noDEG"
      ),
      contrast = !!name
    ) %>% 
    # dplyr::filter(diff != "noDEG") %>% 
    dplyr::select(geneId, !!lfcCol, padj, diff, contrast, !!!otherCols)
  
  return(degs)
}


i <- 1
degData <- NULL

for(i in 1:nrow(rnaseqInfo)){
  dt <- get_foldchange(degFile = rnaseqInfo$deg[i], name = rnaseqInfo$comparison[i],
                       lfcCol = col_lfc, otherCols = c("ENSEMBL", "GENE_NAME", "DESCRIPTION"))
  degData <- dplyr::bind_rows(degData, dt)
}

# geneInfo <- dplyr::left_join(geneInfo, dt, by = c("geneId" = "geneId"))

degData <- tidyr::unite(data = degData, col = "group", contrast, diff, sep = ".", remove = FALSE)


###########################################################################

degLists <- dplyr::filter(degData, diff != "noDEG") %>% 
  dplyr::group_by(contrast, diff, group) %>% 
  dplyr::summarise(
    n = n(),
    geneList = list(geneId))

geneList <- purrr::set_names(x = degLists$geneList, nm = degLists$group)

title <- paste("DEG overlap in", paste(degResults, collapse = ", "), "\n", collapse = " ")
wrap_80 <- scales::wrap_format(80)
title <- wrap_80(title)

cm <- make_comb_mat(geneList, mode = "distinct")

set_name(cm)
set_size(cm)
comb_name(cm)
comb_size(cm)
comb_degree(cm)

grpsDD <- degLists$group[which(degLists$diff == "down")]
combDD <- paste(as.numeric(set_name(cm) %in% grpsDD), collapse = "")
grpsUU <- degLists$group[which(degLists$diff == "up")]
combUU <- paste(as.numeric(set_name(cm) %in% grpsUU), collapse = "")

colorComb <- structure(rep("black", times = length(comb_name(cm))), names = comb_name(cm))
colorComb[c(combDD, combUU)] <- c("red", "blue")

pt <- UpSet(
  m = cm,
  pt_size = unit(7, "mm"), lwd = 3,
  set_order = names(geneList),
  comb_col = colorComb,
  top_annotation = HeatmapAnnotation(
    foo = anno_empty(border = FALSE),
    "combSize" = anno_text(
      x = paste("(", comb_size(cm), ")", sep = ""),
      just = "center", rot = 0
    ),
    "Intersection\nsize" = anno_barplot(
      x = comb_size(cm), 
      border = FALSE, 
      gp = gpar(fill = colorComb), 
      height = unit(2, "cm")
    ),
    annotation_name_side = "left", 
    annotation_name_rot = 0
  ),
  right_annotation = upset_right_annotation(
    m = cm, bar_width = 0.5
  ),
  row_names_max_width = max_text_width(
    set_name(cm), gp = gpar(fontsize = 12)
  ),
  width = unit(12, "cm"), height = unit(6, "cm")
)


pdf(file = paste(outPrefix, ".overlap_upset.pdf", sep = ""), width = 12, height = 8)
# png(filename = paste(outPrefix, ".overlap_upset.png", sep = ""), height = 1500, width = 4000, res = 200)
draw(
  pt,
  column_title = title,
)
dev.off()


###########################################################################
## GO enrichment for individual groups
degGroupsDf <- degData %>%
  pivot_wider(
    id_cols = c(geneId, ENSEMBL, GENE_NAME, DESCRIPTION),
    names_from = c(contrast),
    values_from = c(diff, log2FoldChange, padj),
    names_sep = "."
  )

## save the data
readr::write_tsv(x = degGroupsDf, path = paste(outPrefix, ".DEG_overlap.tab", sep = ""))


## remove the genes which are "noDEG" in both conditions to perform GO enrichment
degGroupsDf <- dplyr::filter_at(
  .tbl = degGroupsDf,
  .vars = vars(starts_with("diff")),
  .vars_predicate = any_vars(. != "noDEG")
)


degGroupsGo <- dplyr::group_by_at(degGroupsDf, .vars = vars(starts_with("diff."))) %>%
  # dplyr::summarise(n = n()) %>%
  dplyr::do(
    topGO_enrichment(
      goMapFile = file_topGO,
      genes = .$geneId,
      type = "BP", goNodeSize = 5,
      orgdb = orgDb, keytype = col_degOrgdbKey,
      topgoColumn = col_topGO, geneNameColumn = col_geneName
    )
  )


readr::write_tsv(x = degGroupsGo, path = paste(outPrefix, ".DEG_groups.topGO.tab", sep = ""))

###########################################################################
## plot GO enrichment figures
degGroupsGo <- suppressMessages(
  readr::read_tsv(file = paste(outPrefix, ".DEG_groups.topGO.tab", sep = ""))
)


goBarPlots <- degGroupsGo %>%
  dplyr::group_by_at(.vars = vars(starts_with("diff."))) %>%
  dplyr::slice(1:10) %>%
  dplyr::mutate(
    groupId = dplyr::cur_group_id()
  ) %>%
  tidyr::unite(col = "group", groupId, starts_with("diff."), remove = FALSE) %>%
  dplyr::do(
    plots = enrichment_bar(
      df = ., title = paste(title, ":", unique(.$group)),
      pvalCol = "weightedFisher",
      termCol = "Term",
      colorCol = "group",
      countCol = "Significant"
    ),
    group = unique(.$group)
  ) %>%
  tidyr::unnest(group)


goPlotList <- goBarPlots$plots %>%
  purrr::set_names(goBarPlots$group)



goPlotList <- cowplot::align_plots(plotlist = goPlotList, align = "v", axis = "rl")

pdf(file = paste(outPrefix, ".DEG_overlap.topGO.barplot.pdf", sep = ""),
    width = 10, height = 4, onefile = TRUE)

for (i in names(goPlotList)) {
  plot(goPlotList[[i]])
}

dev.off()


###########################################################################




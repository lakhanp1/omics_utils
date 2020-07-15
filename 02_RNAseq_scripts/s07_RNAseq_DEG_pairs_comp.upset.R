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
## This script plots the combined upset plots for pairwise DEG comparisons
##

rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")
source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/s01_topGO_functions.R")

###########################################################################
analysisName <- "DEG_pair_comparison"

outDir <- here::here("analysis", "04_DEG_compare")
outPrefix <- paste(outDir, analysisName, sep = "/")


file_RNAseq_info <- here::here("data", "reference_data", "DESeq2_DEG_info.txt")
diffDataPath <- here::here("analysis", "02_DESeq2_diff")

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

file_config <- here::here("analysis", "04_DEG_compare", "DEG_pair_compare.conf.tab")

###########################################################################

paircompConf <- suppressMessages(readr::read_tsv(file = file_config))

degResults <-  union(paircompConf$deg1, paircompConf$deg2)

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degResults)

rnaseqInfoList <- purrr::transpose(rnaseqInfo) %>% 
  purrr::set_names(nm = map(., "comparison"))

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


degLists <- purrr::map(
  .x = rnaseqInfoList,
  .f = function(x){
    
    dt <- get_foldchange(degFile = x$deg, name = x$comparison,
                         lfcCol = col_lfc, otherCols = c("ENSEMBL")) %>% 
      dplyr::filter(diff != "noDEG") %>% 
      dplyr::mutate(
        drug = x$drug,
        concentration = x$concentration,
        time = x$time
      ) %>% 
      tidyr::unite(col = "group", sep = ".", drug, diff, remove = FALSE)

    split(x = dt$geneId, f = dt$group)
    
  }
)


cmList <- purrr::transpose(paircompConf) %>%
  purrr::set_names(nm = map(., "degPairId")) %>% 
  purrr::map(
    .f = function(x){
      cm <- make_comb_mat(c(degLists[[x$deg1]], degLists[[x$deg2]]), mode = "distinct")
    }
  )

sapply(cmList, comb_size)
cmNormList <- normalize_comb_mat(cmList)
sapply(cmNormList, comb_size)
sapply(cmNormList, set_name)
sapply(cmNormList, set_size)
sapply(cmNormList, comb_name)
sapply(cmNormList, comb_degree)


tmpCm <- cmNormList[[1]]
set_name(tmpCm)
set_size(tmpCm)
comb_name(tmpCm)
comb_size(tmpCm)
comb_degree(tmpCm)

## identify the up-up and down-down combinations and use different color
grpsDD <- grepl(pattern = ".down", set_name(tmpCm))
combDD <- paste(as.numeric(grpsDD), collapse = "")
grpsUU <- grepl(pattern = ".up", set_name(tmpCm))
combUU <- paste(as.numeric(grpsUU), collapse = "")

colorComb <- structure(rep("black", times = length(comb_name(tmpCm))), names = comb_name(tmpCm))
colorComb[c(combDD, combUU)] <- c("blue", "red")

## show the up-up and down-down combination first 
combOrder <- forcats::as_factor(comb_name(tmpCm)) %>% 
  forcats::fct_relevel(combDD, combUU) %>% 
  order()


ht_list <- NULL


for (i in seq_along(cmNormList)) {
  
  ## generate Upset plot
  pt <- UpSet(
    m = cmNormList[[i]],
    pt_size = unit(7, "mm"), lwd = 3,
    set_order = set_name(tmpCm),
    comb_order = combOrder,
    comb_col = colorComb,
    row_title = names(cmNormList)[i],
    top_annotation = HeatmapAnnotation(
      foo = anno_empty(border = FALSE),
      "combSize" = anno_text(
        x = paste("(", comb_size(cmNormList[[i]]), ")", sep = ""),
        just = "center", rot = 0
      ),
      "Intersection\nsize" = anno_barplot(
        x = comb_size(cmNormList[[i]]), 
        border = FALSE, 
        gp = gpar(fill = colorComb), 
        height = unit(4, "cm")
      ),
      annotation_name_side = "left", 
      annotation_name_rot = 0
    ),
    right_annotation = NULL,
    # right_annotation = upset_right_annotation(
    #   m = cmNormList[[i]], bar_width = 0.5
    # ),
    row_names_max_width = max_text_width(
      set_name(cmNormList[[i]]), gp = gpar(fontsize = 12)
    ),
    width = unit(12, "cm"), height = unit(3, "cm")
  )
  
  ht_list <- ht_list %v% pt
}


pdf(file = paste(outPrefix, ".combined_upset.pdf", sep = ""), width = 8, height = 14)
# png(filename = paste(outPrefix, ".overlap_upset.png", sep = ""), height = 1500, width = 4000, res = 200)
draw(ht_list)
dev.off()















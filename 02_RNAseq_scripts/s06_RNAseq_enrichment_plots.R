suppressPackageStartupMessages(library(tidyverse))

## This script
## 1) read the tabular config file for different genesets
##
## Config file format: TAB delited file
## <deg>  <termId>  <enrichmentType>  <degCategory>  <title>  <output>
## deg: DEG ID from RNAseq_info
## termId: a ; separated GO term or KEGG pathway IDs to be plotted
## enrichmentType: one of topGO/keggProfile
## degCategory: one of up/down
## title: plot title
## output: a suffix for output file
##
## 2) generate the functional enrichment scatter plots


rm(list = ls())
source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/s01_topGO_functions.R")
source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

####################################################################

diffDataPath <- here::here("analysis", "02_DESeq2_diff")
file_termSet <- here::here("analysis", "02_DESeq2_diff", "enrichment_scatter.config.tab")

file_RNAseq_info <- here::here("data", "reference_data", "DESeq2_DEG_info.txt")
file_sampleInfo <- here::here("data", "reference_data", "sample_info.txt")


# analysisName <- "enrichment_plots"
# degResult <- "DKO_vs_WT"
# 
# diffDataPath <- here::here("analysis", "02_DESeq2_diff")
# outDir <- here::here("analysis", "02_DESeq2_diff", degResult, analysisName)
# 
# if(!dir.exists(outDir)){
#   dir.create(path = outDir)
# }
# 
# outPrefix <- paste(outDir, "/", degResult, ".enrichment_scatter", sep = "")
# 

###########################################################################

colList <- list(
  topGO = list(
    col_pval = "pvalue", col_richness = "richness",
    col_geneCount = "Significant",
    col_term = "Term", col_id = "GO.ID",
    title = "topGO enrichment"
  ),
  keggProfile = list(
    col_pval = "pvalue", col_richness = "richness",
    col_geneCount = "Significant",
    col_term = "Pathway_Name", col_id = "pathway_id",
    title = "keggProfile enrichment"
  )
)


termSet <- suppressMessages(readr::read_tsv(file = file_termSet))

termSet <- dplyr::mutate(
  termSet,
  termId = stringr::str_split(string = termId, pattern = ";")
) %>% 
  dplyr::mutate(
    termId = purrr::map(termId, unique)
  )

# rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
#   dplyr::filter(comparison == degResult)
# 

###########################################################################

setRow <- 2

degResult <- termSet$deg[setRow]
outDir <- paste(diffDataPath, "/", degResult, "/geneset_plots", sep = "")

if(!dir.exists(outDir)){
  dir.create(path = outDir)
}

outPrefix <- paste(outDir, "/", degResult, ".enrichment_scatter.", sep = "")


rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == termSet$deg[setRow])


enrichmentType <- termSet$enrichmentType[setRow]
degCategories <- unlist(stringr::str_split(string = termSet$degCategory[setRow], pattern = ";"))
plotTitle <- termSet$title[setRow]
plotOutSuffix <- termSet$output[setRow]


enrichmentData <- suppressMessages(readr::read_tsv(rnaseqInfo[[enrichmentType]])) %>% 
  dplyr::mutate(comparison = rnaseqInfo$comparison)


enrichmentSet <- tibble::tibble(termId = unlist(termSet$termId[setRow])) %>% 
  dplyr::left_join(
    y = enrichmentData,
    by = structure(c(colList[[enrichmentType]]$col_id), names = c("termId"))
  ) %>% 
  dplyr::filter(category %in% degCategories)


pt <- enrichment_scatter(
  df = enrichmentSet,
  title = paste(colList[[enrichmentType]]$title, ":", plotTitle),
  pvalCol = colList[[enrichmentType]]$col_pval,
  termCol = colList[[enrichmentType]]$col_term,
  xVar = colList[[enrichmentType]]$col_richness,
  sizeVar = colList[[enrichmentType]]$col_geneCount
)


# ## optional facetting for multiple RNAseq enrichment results
# pt <- pt +
#   facet_grid(facets = category ~ ., scales = "free_y") +
#   theme(strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 14, face = "bold"))


ggsave(
  filename = paste(outPrefix, ".", plotOutSuffix, ".pdf", sep = ""),
  plot = pt,
  width = 11, height = 6, units = "in"
)



library(tidyverse)


## generate enrichment plots for available functional enrichment results


rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")

###########################################################################


analysisName <- "PG_HE_8mpf_vs_PG_WT_8mpf"
enrichmentType <- "keggProfile"

degResults <- c(analysisName)

outDir <- here::here("analysis", "02_DESeq2_diff", analysisName)
outPrefix <- paste(outDir, "/", analysisName, ".", enrichmentType, sep = "")
file_RNAseq_info <- here::here("data", "RNAseq_info.txt", sep = "")

file_termSet <- here::here("analysis", "02_DESeq2_diff", "combined_summary", "pathway_interest.txt")


colList <- list(
  topGO = list(
    col_pval = "weightedFisher", col_richness = "richness",
    col_term = "Term", col_id = "GO.ID",
    title = "topGO enrichment"
  ),
  keggProfile = list(
    col_pval = "pvalue", col_richness = "richness",
    col_term = "Pathway_Name", col_id = "pathway_id",
    title = "keggProfile enrichment"
  )
)


termSet <- suppressMessages(readr::read_tsv(file = file_termSet))

###########################################################################

rnaseqInfo <- suppressMessages(readr::read_tsv(file = file_RNAseq_info)) %>% 
  dplyr::filter(comparison %in% degResults)


enrichmentData <- NULL

for (i in 1:nrow(rnaseqInfo)) {
  df <- suppressMessages(readr::read_tsv(rnaseqInfo[[enrichmentType]][i])) %>% 
    dplyr::mutate(comparison = rnaseqInfo$comparison[i])
  
  enrichmentData <- dplyr::bind_rows(enrichmentData, df)
}

enrichmentData <- dplyr::filter(
  enrichmentData, !!sym(colList[[enrichmentType]]$col_id) %in% !!termSet$id
)

pt <- enrichment_scatter(df = enrichmentData,
                         title = paste(colList[[enrichmentType]]$title, analysisName),
                         pvalCol = colList[[enrichmentType]]$col_pval,
                         termCol = colList[[enrichmentType]]$col_term,
                         xVar = "category",
                         sizeVar = colList[[enrichmentType]]$col_richness)


# pt <- pt +
#   facet_grid(facets = . ~ comparison) +
#   theme(strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 14, face = "bold"))

png(filename = paste(outPrefix, ".png", sep = ""),
    width = 2500, height = 2500, res = 250)

pt
dev.off()

pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 12, height = 8)
pt
dev.off()







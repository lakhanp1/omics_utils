library(tidyverse)
library(data.table)
library(org.HSapiens.gencodev30.eg.db)
library(fgsea)
library(msigdbr)
library(DT)
library(clusterProfiler)

## RNAseq DEG functional analysis using clusterProfiler

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")

###########################################################################

analysisName <- "MHCC97L_AB_vs_B"
outDir <- here::here("analysis", "02_DESeq2_diff", analysisName)
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_diff <- paste(outDir, "/", analysisName, ".DEG_all.txt", sep = "")

orgDb <- org.HSapiens.gencodev30.eg.db

cutoff_qval <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

# m_df <- msigdbr(species = "Homo sapiens")
# m_list <- dplyr::filter(m_df, gs_cat %in% c("H", "C2")) %>% 
#   split(x = .$gene_symbol, f = .$gs_name)

entrezGids <- AnnotationDbi::select(x = orgDb,
                                    columns = c("NCBI_ID"),
                                    keys = keys(orgDb)) %>% 
  dplyr::filter(!is.na(NCBI_ID))

###########################################################################

degs <- suppressMessages(readr::read_tsv(file = file_diff)) %>% 
  dplyr::mutate(rankMetric = (-log10(padj) * shrinkLog2FC)) %>% 
  dplyr::filter(!is.na(rankMetric)) %>% 
  dplyr::arrange(desc(rankMetric)) %>% 
  dplyr::left_join(y = entrezGids, by = c("ENSEMBL" = "GID"))

geneList <- tibble::deframe(x = dplyr::select(degs, GENE_NAME, rankMetric))

downDegs <- dplyr::filter(degs, padj <= cutoff_qval & shrinkLog2FC <= cutoff_down)
upDegs <- dplyr::filter(degs, padj <= cutoff_qval & shrinkLog2FC >= cutoff_up)

contrast <- unique(upDegs$contrast)
###########################################################################
## clusterProfiler: GO enrichment
ego_up <- enrichGO(gene = unique(upDegs$ENSEMBL),
                   universe = na.omit(unique(degs$ENSEMBL)),
                   OrgDb = orgDb,
                   ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   keyType = "ENSEMBL",
                   readable = FALSE)

ego_up <- simplify(x = ego_up)

ego_down <- enrichGO(gene = unique(downDegs$ENSEMBL),
                     universe = na.omit(unique(degs$ENSEMBL)),
                     OrgDb = orgDb,
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     keyType = "ENSEMBL",
                     readable = FALSE)

ego_down <- simplify(x = ego_down)

ego <- dplyr::bind_rows(
  dplyr::mutate(.data = as.data.frame(ego_up), category = "up"),
  dplyr::mutate(.data = as.data.frame(ego_down), category = "down")
) %>% 
  dplyr::mutate(contrast = contrast)


readr::write_tsv(x = ego, path = paste(outPrefix, ".clusterProfiler.GO.tab", sep = ""))

###########################################################################
## clusterProfiler: KEGG pathway enrichment
ekegg_up <- enrichKEGG(gene = na.omit(unique(upDegs$NCBI_ID)),
                       organism = 'hsa',
                       pvalueCutoff = 0.05)

ekegg_down <- enrichKEGG(gene = na.omit(unique(downDegs$NCBI_ID)),
                         organism = 'hsa',
                         pvalueCutoff = 0.05)

ekegg <- dplyr::bind_rows(
  dplyr::mutate(.data = as.data.frame(ekegg_up), category = "up"),
  dplyr::mutate(.data = as.data.frame(ekegg_down), category = "down")
) %>% 
  dplyr::mutate(contrast = contrast)


readr::write_tsv(x = ekegg, path = paste(outPrefix, ".clusterProfiler.kegg.tab", sep = ""))


## up and down DEG
cp_kegg <- compareCluster(
  geneClusters = list(up = na.omit(unique(upDegs$NCBI_ID)),
                      down = na.omit(unique(downDegs$NCBI_ID))),
  fun = "enrichKEGG",
  organism     = 'hsa'
)

gg_cp_kegg <- dotplot(cp_kegg, showCategory = 100) +
  labs(title = paste(analysisName, "KEGG pathway enrichment")) +
  theme(
    plot.title = element_text(hjust = 1)
  )


png(filename = paste(outPrefix, ".clusterProfiler.kegg.png", sep = ""),
    width = 1500, height = 1500, res = 200)

gg_cp_kegg

dev.off()







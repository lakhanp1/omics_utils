library(tidyverse)
library(data.table)
library(fgsea)
library(msigdbr)
library(DT)
library(clusterProfiler)
library(grid)
require(XLConnect)
options(java.parameters = "- Xmx4g")
xlcFreeMemory()
library(org.AFumigatus293.eg.db)


## RNAseq DEG functional analysis using clusterProfiler

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/topGO_functions.R")
source("E:/Chris_UM/GitHub/omics_util/02_RNAseq_scripts/s02_DESeq2_functions.R")

###########################################################################

degResult <- "5A9_AA_vs_5A9_C"

file_RNAseq_info <- here::here("data", "RNAseq_info.txt")
diffDataPath <- here::here("analysis", "01_RNAseq_data")

outDir <- here::here("analysis", "01_RNAseq_data", degResult)
outPrefix <- paste(outDir, "/", degResult, sep = "")

# file_msigDesc <- "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/msigDB_geneset_desc.tab"

orgDb <- org.AFumigatus293.eg.db
keggOrg <- 'afm'
keggIdCol <- "KEGG_ID"
file_topGO <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/annotation_resources/geneid2go.AFumigatus_Af293.topGO.map"

cutoff_qval <- 0.05
cutoff_lfc <- 0.585
cutoff_up <- cutoff_lfc
cutoff_down <- -1 * cutoff_lfc

###########################################################################

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison == degResult)

degs <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[1])) %>% 
  dplyr::mutate(rankMetric = (-log10(pvalue) * sign(shrinkLog2FC))) %>%
  dplyr::arrange(desc(rankMetric)) %>% 
  dplyr::filter(!is.na(rankMetric))

if(! keggIdCol %in% colnames(degs)){
  keggInfo <- suppressMessages(
    AnnotationDbi::select(x = orgDb, keys = degs$geneId,
                                    keytype = "GID", columns = keggIdCol)
  ) %>% 
    dplyr::filter(!is.na(!!sym(keggIdCol))) %>% 
    dplyr::rename(geneId = GID)
  
  degs <- dplyr::left_join(x = degs, y = keggInfo, by = "geneId")
}
  
downDegs <- dplyr::filter(degs, padj <= cutoff_qval & shrinkLog2FC <= cutoff_down)
upDegs <- dplyr::filter(degs, padj <= cutoff_qval & shrinkLog2FC >= cutoff_up)

contrast <- unique(upDegs$contrast)

geneList <- dplyr::filter(degs, !is.na(geneId)) %>% 
  dplyr::select(geneId, rankMetric) %>% 
  tibble::deframe()

###########################################################################
# ## clusterProfiler: GO enrichment
# ego_up <- enrichGO(gene = unique(upDegs$geneId),
#                    OrgDb = orgDb,
#                    ont = "BP", pAdjustMethod = "BH",
#                    pvalueCutoff = 0.05,
#                    qvalueCutoff = 0.05,
#                    keyType = "ENSEMBL",
#                    readable = FALSE)
# 
# ego_up <- simplify(x = ego_up)
# 
# ego_down <- enrichGO(gene = unique(downDegs$geneId),
#                      OrgDb = orgDb,
#                      ont = "BP", pAdjustMethod = "BH",
#                      pvalueCutoff = 0.05,
#                      qvalueCutoff = 0.05,
#                      keyType = "ENSEMBL",
#                      readable = FALSE)
# 
# ego_down <- simplify(x = ego_down)
# 
# ego_res <- dplyr::bind_rows(
#   dplyr::mutate(.data = as.data.frame(ego_up), category = "up"),
#   dplyr::mutate(.data = as.data.frame(ego_down), category = "down")
# ) %>%
#   dplyr::mutate(contrast = contrast)
# 
# 
# readr::write_tsv(x = ego_res, path = paste(outPrefix, ".clusterProfiler.GO.tab", sep = ""))


## topGO GO enrichment
topgo_up <- topGO_enrichment(goMapFile = file_topGO,
                             genes = unique(upDegs$geneId),
                             type = "BP",
                             goNodeSize = 5)

topgo_down <- topGO_enrichment(goMapFile = file_topGO,
                               genes = unique(downDegs$geneId),
                               type = "BP",
                               goNodeSize = 5)

topgo_res <- dplyr::bind_rows(
  dplyr::mutate(.data = as.data.frame(topgo_up), category = "up"),
  dplyr::mutate(.data = as.data.frame(topgo_down), category = "down")
) %>%
  dplyr::mutate(contrast = contrast)

readr::write_tsv(x = topgo_res, path = paste(outPrefix, ".topGO.tab", sep = ""))

###########################################################################
# ## clusterProfiler: KEGG pathway enrichment
# ekegg_up <- enrichKEGG(gene = na.omit(unique(upDegs$NCBI_ID)),
#                        organism = keggOrg,
#                        pvalueCutoff = 0.05)
# 
# ekegg_down <- enrichKEGG(gene = na.omit(unique(downDegs$NCBI_ID)),
#                          organism = keggOrg,
#                          pvalueCutoff = 0.05)
# 
# ekegg_res <- dplyr::bind_rows(
#   dplyr::mutate(.data = as.data.frame(ekegg_up), category = "up"),
#   dplyr::mutate(.data = as.data.frame(ekegg_down), category = "down")
# ) %>%
#   dplyr::mutate(contrast = contrast)
# 
# 
# readr::write_tsv(x = ekegg_res, path = paste(outPrefix, ".clusterProfiler.kegg.tab", sep = ""))


# ## up and down DEG
# cp_kegg <- compareCluster(
#   geneClusters = list(up = na.omit(unique(upDegs$NCBI_ID)),
#                       down = na.omit(unique(downDegs$NCBI_ID))),
#   fun = "enrichKEGG",
#   organism = keggOrg
# )
# 
# gg_cp_kegg <- dotplot(cp_kegg, showCategory = 100) +
#   labs(title = paste(analysisName, "KEGG pathway enrichment")) +
#   theme(
#     plot.title = element_text(hjust = 1)
#   )
# 
# 
# png(filename = paste(outPrefix, ".clusterProfiler.kegg.png", sep = ""),
#     width = 1500, height = 1500, res = 200)
# 
# gg_cp_kegg
# 
# dev.off()


## KEGGprofile::find_enriched_pathway
keggp_up <- keggprofile_enrichment(
  genes = as.character(na.omit(unique(upDegs[[keggIdCol]]))), orgdb = orgDb,
  keytype = keggIdCol, keggIdCol = keggIdCol, keggOrg = keggOrg
)

keggp_down <- keggprofile_enrichment(
  genes = as.character(na.omit(unique(downDegs[[keggIdCol]]))), orgdb = orgDb,
  keytype = keggIdCol, keggIdCol = keggIdCol, keggOrg = keggOrg
)

keggp_res <- dplyr::bind_rows(
  dplyr::mutate(.data = as.data.frame(keggp_up), category = "up"),
  dplyr::mutate(.data = as.data.frame(keggp_down), category = "down")
) %>%
  dplyr::mutate(contrast = contrast)

readr::write_tsv(x = keggp_res, path = paste(outPrefix, ".keggProfile.tab", sep = ""))

# ###########################################################################
# ## GSEA
# msigdbr_show_species()
# msig_df <- msigdbr(species = "Danio rerio", category = c("H", "C2", "C5")) %>% 
#   dplyr::filter(! gs_subcat %in% c("MF", "CC"))
# 
# msig_list <- split(x = msig_df$entrez_gene, f = msig_df$gs_name)
# 
# length(intersect(names(geneList), unique(msig_df$entrez_gene)))
# 
# vn <- VennDiagram::venn.diagram(
#   x = list(geneList = names(geneList), msig = unique(msig_df$entrez_gene)),
#   filename = NULL,
#   print.mode = c("raw", "percent"),
#   scaled = FALSE
# )
# dev.off()
# grid.draw(vn)
# 
# # msigDescDf <- suppressMessages(readr::read_tsv(file = file_msigDesc))
# # msigDesc <- split(x = msigDescDf$DESCRIPTION_BRIEF, f = msigDescDf$STANDARD_NAME)
# 
# egsea <- GSEA(geneList = geneList,
#               nPerm = 10000,
#               pvalueCutoff = 0.1,
#               minGSSize = 10, maxGSSize = Inf,
#               TERM2GENE = dplyr::select(msig_df, gs_name, entrez_gene))
# 
# egseaDf <- as_tibble(egsea) %>% 
#   dplyr::left_join(y = msigDescDf, by = c("ID" = "STANDARD_NAME")) %>% 
#   dplyr::mutate(contrast = contrast) %>% 
#   dplyr::select(ID, contrast, everything(), -Description)
# 
# readr::write_tsv(x = egseaDf,
#                  path = paste(outPrefix, ".clusterProfiler.GSEA.tab", sep = ""))
# 
# 
# 
# 
# ## plotting specific genesets
# genesetSub <- c("GO_CELL_CYCLE",
#                 "GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
#                 "GO_DNA_REPLICATION")
# 
# plotList <- list()
# 
# pdf(file = paste(outPrefix, ".clusterProfiler.GSEA_enrichmentPlot.pdf", sep = ""),
#     width = 10, height = 8, onefile = TRUE)
# 
# for (setId in genesetSub) {
#   
#   if(setId %in% egsea$ID){
#     pt <- enrichplot::gseaplot2(egsea, geneSetID = setId)
#     wrap_100 <- wrap_format(120)
#     
#     plotSubTitle <- paste(
#       "p-value = ", sprintf(fmt = "%.2E", egseaDf$pvalue[which(egseaDf$ID == setId)]),
#       "; q-value = ", sprintf(fmt = "%.2E", egseaDf$qvalues[which(egseaDf$ID == setId)]),
#       "\n", wrap_100(x = msigDesc[[setId]]),
#       sep = "")
#     
#     pt <- pt +
#       labs(
#         title = paste(setId, ": ", contrast, sep = ""),
#         subtitle = plotSubTitle
#       ) +
#       theme_bw() +
#       theme(panel.grid = element_blank(),
#             panel.border = element_blank())
#     
#     plotList[[setId]] <- pt
#     
#     plot(pt)
#     
#   }
#   
# }
# 
# dev.off()
# 
# 
# ###########################################################################
# ## GSEA enrichment using fgsea
# gseaRes <- fgsea(pathways = msig_list, stats = geneList, nperm = 10000)
# 
# pt2 <- plotEnrichment(pathway = msig_list[[setId]],
#                       stats = geneList) +
#   labs(
#     title = paste(setId, ":", contrast),
#     subtitle = wrap_100(x = msigDesc[[setId]]),
#     x = "Rank in ordered dataset",
#     y = "Enrichment Score") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14, face = "bold"))
# 
# 


# ###########################################################################
excelOut <- paste(outPrefix, ".enrichment.xlsx", sep = "")
unlink(excelOut, recursive = FALSE, force = FALSE)
exc = loadWorkbook(excelOut , create = TRUE)
xlcFreeMemory()

wrkSheet <- "topGO"
createSheet(exc, name = wrkSheet)
createFreezePane(exc, sheet = wrkSheet, 2, 2)
writeWorksheet(object = exc, data = topgo_res, sheet = wrkSheet)
setAutoFilter(object = exc, sheet = wrkSheet,
              reference = aref(topLeft = "A1", dimension = dim(topgo_res)))

wrkSheet <- "keggProfile"
createSheet(exc, name = wrkSheet)
createFreezePane(exc, sheet = wrkSheet, 2, 2)
writeWorksheet(object = exc, data = keggp_res, sheet = wrkSheet)
setAutoFilter(object = exc, sheet = wrkSheet,
              reference = aref(topLeft = "A1", dimension = dim(keggp_res)))

# wrkSheet <- "clusterProfiler_GO"
# createSheet(exc, name = wrkSheet)
# createFreezePane(exc, sheet = wrkSheet, 2, 2)
# writeWorksheet(object = exc, data = ego_res, sheet = wrkSheet)
# setAutoFilter(object = exc, sheet = wrkSheet,
#               reference = aref(topLeft = "A1", dimension = dim(ego_res)))
# 
# wrkSheet <- "clusterProfiler_KEGG"
# createSheet(exc, name = wrkSheet)
# createFreezePane(exc, sheet = wrkSheet, 2, 2)
# writeWorksheet(object = exc, data = ekegg_res, sheet = wrkSheet)
# setAutoFilter(object = exc, sheet = wrkSheet,
#               reference = aref(topLeft = "A1", dimension = dim(ekegg_res)))

setColumnWidth(object = exc, sheet = 1:2, column = 1, width = -1)
setColumnWidth(object = exc, sheet = 1:2, column = 2, width = c(13000))

xlcFreeMemory()
saveWorkbook(exc)





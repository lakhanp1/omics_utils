library(chipmine)
library(org.Anidulans.eg.db)

## this script generates profile plot for single sample

rm(list = ls())

path = "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/TF_data/An_kdmB_20h_HA_1"
setwd(path)


## genes to read
file_exptInfo <-"E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/sampleInfo.txt"
file_genes <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/referenceData/AN_genesForPolII.bed"
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"


orgDb <- org.Anidulans.eg.db

outpath <- path
sampleId <- "An_kdmB_20h_HA_1"
file_mat <- "An_kdmB_20h_HA_1_normalized_profile.tab.gz"
file_bwMat <- "An_kdmB_20h_HA_1_normalizedMatrix.tab.gz"
file_bw <- "An_kdmB_20h_HA_1_normalized.bw"
outPrefix <- "test_An_kdmB_20h_HA_1"
file_cluster <- "An_kdmB_20h_HA_1.kmeans.clusters.txt"
doClustering <- FALSE

##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

## gene information annotations: cluster and TF and polII expression values
geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = geneSet)

head(geneInfo)

##################################################################################


## read the profile matrix
mat1 <- import_profile_from_file(file = file_mat,
                                 source = "deeptools",
                                 signalName = sampleId,
                                 selectGenes = geneInfo$gene)


mat2 <- import_profile_from_file(file = file_bwMat,
                                 source = "normalizedmatrix",
                                 signalName = sampleId,
                                 selectGenes = geneInfo$gene)


## check the distribution in data
quantile(mat1, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
col_fun <- colorRamp2(quantile(mat1, c(0.50, 0.995), na.rm = T), c("white", "red"))

## get first primary profile heatmap
profile1 <- profile_heatmap(profileMat = mat1,
                            signalName = sampleId,
                            profileColor = col_fun,
                            geneGroups = file_cluster,
                            createClusters = doClustering,
                            ylimFraction = 0.995)


anGl <- gene_length_heatmap_annotation(bedFile = file_genes, genes = geneInfo$gene)


htlist <- profile1$rowAnno + anGl$an + profile1$heatmap

# draw Heatmap and add the annotation name decoration
png(filename = paste0(outPrefix, ".png", collapse = ""), width=4000, height=6000, res = 450)

draw(htlist,
     main_heatmap = sampleId,
     # annotation_legend_list = list(profile1$legend),
     column_title = sampleId,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     # gap = unit(c(2), "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()


##################################################################################

newDf <- dplyr::filter(profile1$cluster, cluster %in% c("Cluster_1", "Cluster_4", "Cluster_3")) %>% 
  dplyr::select(gene, cluster)


## read the profile matrix
subMat <- mat1[newDf$gene, ]


subProfile <- profile_heatmap(profileMat = subMat,
                              signalName = sampleId,
                              profileColor = col_fun,
                              geneGroups = newDf,
                              ylimFraction = profile1$ylim)


subGl <- gene_length_heatmap_annotation(bedFile = file_genes, genes = newDf$gene)


subHtList <- subProfile$rowAnno + subGl$an + subProfile$heatmap

# draw Heatmap and add the annotation name decoration
png(filename = paste0(outPrefix, "_sub.png", collapse = ""), width=4000, height=6000, res = 450)

draw(subHtList,
     main_heatmap = sampleId,
     # annotation_legend_list = list(profile1$legend),
     column_title = "Profile heatmap with subset of gene",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     # gap = unit(c(2), "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()





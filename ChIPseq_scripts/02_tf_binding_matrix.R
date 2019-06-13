library(chipmine)
library(here)
library(FactoMineR)
library(factoextra)
library(NbClust)

rm(list = ls())

outDir <- here::here("analysis", "peak_occupancy")


if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")

TF_dataPath <- here::here("data", "ChIPseq_data")

outPrefix <- paste(outDir, "/ctcf_ChIP", sep = "")

##################################################################################
## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   dataPath = TF_dataPath)

exptDataList <- purrr::transpose(exptData) %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

##################################################################################
## merge the peaks and generate the combination matrix for presence of peak
mat <- combinatorial_binding_matrix(sampleInfo = exptData)

readr::write_tsv(x = mat, path = paste(outPrefix, ".binding_matrix.tab", sep = ""))

##################################################################################
## heatmap
peakMat <- dplyr::select(mat, name, starts_with("overlap")) %>% 
  dplyr::mutate_if(is.logical, as.numeric) %>%
  tibble::column_to_rownames(var = "name") %>% 
  as.matrix()

colnames(peakMat) <- gsub(pattern = "overlap.", replacement = "", x = colnames(peakMat))


enrichmentMat <- dplyr::select(mat, name, starts_with("peakEnrichment")) %>% 
  tibble::column_to_rownames(var = "name") %>% 
  as.matrix() %>% 
  scale_matrix_columns(add_attr = F)

colnames(enrichmentMat) <- gsub(pattern = "peakEnrichment.", replacement = "", x = colnames(enrichmentMat))

## transpose the matrices to plot data horizontally
peakMat <- t(peakMat)
enrichmentMat <- t(enrichmentMat)

quantile(enrichmentMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

colorEnrichment <- colorRamp2(
  # breaks = c(0, quantile(enrichmentMat, c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9), na.rm = T)),
  breaks = c(-2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 5, 7),
  colors = c(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn"), "blue", "darkblue")
)


## use top 10000 regions for sample clustering
## row clustering for samples
rowClust <- hclust(dist(peakMat[, sample.int(n = ncol(peakMat), size = 10000)]))
## order of columns: for peaks
colomnOrder <- order(matrixStats::colMeans2(enrichmentMat, na.rm = T))

peakHt <- Heatmap(
  matrix = peakMat,
  name = "peak",
  col = c("1" = "#1f78b4", "0" = "gray95"),
  row_title = "CTCF peak occupancy",
  column_title = NULL,
  column_km = 12,
  cluster_rows = rowClust,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(title = NULL, at = c(1), labels = "CTCF peak  ",
                              direction = "horizontal",
                              grid_height = unit(1,"cm"),
                              grid_width = unit(1, "cm"),
                              labels_gp = gpar(fontsize = 18))
)


enrichmentHt <- Heatmap(
  matrix = enrichmentMat,
  name = "enrichment",
  col = colorEnrichment,
  na_col = "gray95",
  column_title = NULL,
  row_title = "CTCF peak enrichment score",
  cluster_rows = rowClust,
  cluster_columns = FALSE,
  show_column_names = FALSE, 
  heatmap_legend_param = list(title = "z-score(peak enrichment)",
                              direction = "horizontal",
                              title_position = "topcenter",
                              title_gp = gpar(fontsize = 18),
                              legend_width = unit(7, "cm"),
                              labels_gp = gpar(fontsize = 18))
)

htList <- peakHt %v% enrichmentHt

png(filename = paste(outPrefix, ".binding_matrix.k12.png", sep = ""), width = 6000, height = 3500, res = 350)
draw(
  htList,
  column_order = colomnOrder,
  gap = unit(3, "mm"),
  column_title = paste("CTCF peaks\n",
                       "top heatmap: peak occupancy; bottom heatmap: macs2 peakEnrichment"),
  heatmap_legend_side = "bottom"
)
dev.off()


# combMat <- make_comb_mat(peakMat, mode = "distinct")
##################################################################################
## chromosome wise plots

for(chr in unique(mat$seqnames)){
  
  peakMat <- dplyr::filter(mat, seqnames == !!chr) %>% 
    dplyr::arrange(start) %>% 
    dplyr::select(name, starts_with("overlap")) %>% 
    dplyr::mutate_if(is.logical, as.numeric) %>%
    tibble::column_to_rownames(var = "name") %>% 
    as.matrix()
  
  colnames(peakMat) <- gsub(pattern = "overlap.", replacement = "", x = colnames(peakMat))
  
  enrichmentMat <- dplyr::filter(mat, seqnames == !!chr) %>% 
    dplyr::arrange(start) %>% 
    dplyr::select(name, starts_with("peakEnrichment")) %>% 
    tibble::column_to_rownames(var = "name") %>% 
    as.matrix() %>% 
    scale_matrix_columns(add_attr = F)
  
  colnames(enrichmentMat) <- gsub(pattern = "peakEnrichment.", replacement = "", x = colnames(enrichmentMat))
  
  
  ## transpose the matrices to plot data horizontally
  peakMat <- t(peakMat)
  enrichmentMat <- t(enrichmentMat)
  
  quantile(enrichmentMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
  
  colorEnrichment <- colorRamp2(
    # breaks = c(0, quantile(enrichmentMat, c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9), na.rm = T)),
    breaks = c(-2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 5, 7),
    colors = c(RColorBrewer::brewer.pal(n = 11, name = "RdYlGn"), "blue", "darkblue")
  )
  
  
  peakHt <- Heatmap(
    matrix = peakMat,
    name = "peak",
    col = c("1" = "#1f78b4", "0" = "gray95"),
    row_title = "CTCF peak occupancy",
    column_title = NULL,
    cluster_rows = rowClust,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(title = NULL, at = c(1), labels = "CTCF peak  ",
                                direction = "horizontal",
                                grid_height = unit(1,"cm"),
                                grid_width = unit(1, "cm"),
                                labels_gp = gpar(fontsize = 18))
  )
  
  
  enrichmentHt <- Heatmap(
    matrix = enrichmentMat,
    name = "enrichment",
    col = colorEnrichment,
    na_col = "gray95",
    column_title = NULL,
    row_title = "CTCF peak enrichment score",
    cluster_rows = rowClust,
    cluster_columns = FALSE,
    show_column_names = FALSE, 
    heatmap_legend_param = list(title = "z-score(peak enrichment)",
                                direction = "horizontal",
                                title_position = "topcenter",
                                title_gp = gpar(fontsize = 18),
                                legend_width = unit(7, "cm"),
                                labels_gp = gpar(fontsize = 18))
  )
  
  htList <- peakHt %v% enrichmentHt
  
  
  png(filename = paste(outPrefix, ".", chr, ".binding_matrix.png", sep = ""),
      width = 6000, height = 3500, res = 350)
  draw(
    htList,
    gap = unit(3, "mm"),
    column_order = 1:ncol(enrichmentMat),
    column_title = paste(ncol(enrichmentMat), "CTCF peaks on", chr, "(ordered by position)\n",
                         "top heatmap: peak occupancy; bottom heatmap: macs2 peakEnrichment"),
    heatmap_legend_side = "bottom"
  )
  dev.off()
  
}



##################################################################################
## PCA
pcaDf <- dplyr::filter_at(.tbl = mat, .vars = vars(starts_with("overlap.")),
                          .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::select(name, starts_with("peakEnrichment.")) %>%
  data.table::as.data.table() %>%
  melt.data.table(id.vars = "name") %>%
  dcast.data.table(variable ~ name) %>%
  as.data.frame() %>%
  dplyr::rename(sampleId = variable)

pcaDf$sampleId <- gsub(pattern = "peakEnrichment.", replacement = "", x = pcaDf$sampleId)
rownames(pcaDf) <- pcaDf$sampleId

res.pca <- PCA(X = pcaDf, quali.sup = 1)

eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

ind <- get_pca_ind(res.pca)

fviz_pca_ind(res.pca,
             fill.ind = pcaDf$sampleId,
             pointshape = 21,
             repel = TRUE,
             mean.point = FALSE,
             legend.title = "Study",
             pointsize = 3
)







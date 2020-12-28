library(chipmine)
require(XLConnect)
options(java.parameters = "- Xmx4g")
xlcFreeMemory()
library(here)
library(foreach)
library(doParallel)


rm(list = ls())

file_tcp <- here::here("data", "referenceData/tf_ctrl_polII_pairs.txt")

tcpDf <- readr::read_tsv(file = file_tcp, col_names = T, comment = "#")
##################################################################################

cl <- makeCluster(3) #not to overload your computer
registerDoParallel(cl)

i <- 1

foreach(i = 1:nrow(tcpDf),
        .packages = c("chipmine", "XLConnect")) %dopar% {
  
          

  TF_sample <- tcpDf$sampleId[i]
  input_sample <- tcpDf$control[i]
  polII_sample <- tcpDf$polII[i]
  h3_sample <- tcpDf$h3[i]
  
  
  name <- TF_sample
  
  geneFilter <- c("AN5245", "AN3245")
  
  # "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
  matrixType <- "normalizedmatrix_5kb"
  matrixDim = c(500, 200, 100, 10)
  
  tfYlim = 0.996              ##0.999
  # poliiYlim = 0.995           ##0.995
  
  
  file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
  file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
  
  TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
  polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
  hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
  other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")
  
  orgDb <- org.Anidulans.FGSCA4.eg.db
  
  ##################################################################################
  
  outPath <- paste(TF_dataPath, "/", name, sep = "")
  outPrefix_all <- paste(outPath, "/", name, "_allGenes", sep = "")
  outPrefix_expressed <- paste(outPath, "/", name, "_expressedGenes", sep = "")
  outPrefix_sm <- paste(outPath, "/", name, "_SM_genes", sep = "")
  outPrefix_peaks <- paste(outPath, "/", name, "_peaksGenes", sep = "")
  outPrefix_pkExp <- paste(outPath, "/", name, "_pkExpGenes", sep = "")
  
  
  tf_info <- get_sample_information(
    exptInfoFile = file_exptInfo,
    samples = TF_sample,
    dataPath = TF_dataPath,
    profileMatrixSuffix = matrixType
  )
  
  input_info <- get_sample_information(
    exptInfoFile = file_exptInfo,
    samples = input_sample,
    dataPath = TF_dataPath,
    profileMatrixSuffix = matrixType
  )
  
  
  polII_info <- get_sample_information(
    exptInfoFile = file_exptInfo,
    samples = polII_sample,
    dataPath = polII_dataPath,
    profileMatrixSuffix = matrixType
  )
  
  hist_info <- get_sample_information(
    exptInfoFile = file_exptInfo,
    samples = h3_sample,
    dataPath = hist_dataPath,
    profileMatrixSuffix = matrixType
  )
  
  
  sampleInfo <- dplyr::bind_rows(tf_info, input_info, polII_info, hist_info)
  
  
  excelOut <- paste(outPrefix_all, "_clusters", ".xlsx", sep = "")
  
  polIICols <- list(
    exp = structure(polII_sample, names = polII_sample),
    is_expressed = structure(paste("is_expressed", ".", polII_sample, sep = ""), names = polII_sample)
  )
  
  tfCols <- sapply(
    c("peakDist", "featureCovFrac", "hasPeak", "peakCoverage", "peakPosition", "peakId", "peakType",
      "peakPval", "peakEnrichment", "preference", "peakCategory"),
    FUN = function(x){ structure(paste(x, ".", TF_sample, sep = ""), names = TF_sample) },
    simplify = F, USE.NAMES = T)
  
  
  colPal = RColorBrewer::brewer.pal(9, "YlGnBu")
  anWidth = unit(8, "mm")
  
  peakAnArgs <- list()
  anLables <- list()
  
  anLables[[unname(tfCols$peakType)]] <- gsub("peakType", "TSS peak type\n", unname(tfCols$peakType)) %>% gsub("\\(|\\)", "", .)
  anLables[[unname(polIICols$is_expressed)]] <- gsub("is_expressed", "is expressed\n", unname(polIICols$is_expressed)) %>% gsub("\\(|\\)", "", .)
  anLables[["is_SM_gene"]] <- "SM gene"
  anLables[["is_TF"]] <- "Transcription Factor"
  
  profileColors <- list()
  profileYlims <- list()
  
  ##################################################################################
  ## genes to read
  geneSet <-  suppressMessages(
    readr::read_tsv(
      file = file_genes,col_names = c("chr", "start", "end", "geneId", "score", "strand")
    )) %>% 
    dplyr::mutate(length = end - start) %>% 
    dplyr::filter(! geneId %in% geneFilter)
  
  kmClust <- dplyr::left_join(
    x = suppressMessages(readr::read_tsv(file = tf_info$clusterFile[1])),
    y = geneSet, by = "geneId"
  )
  
  ## gene information annotations: cluster and TF and polII expression values
  geneDesc <- AnnotationDbi::select(
    x = orgDb, keys = geneSet$geneId, keytype = "GID",
    columns = c("DESCRIPTION")
  )
  
  smGenes <- AnnotationDbi::select(
    x = orgDb, keys = keys(x = orgDb, keytype = "SM_ID"), keytype = "SM_ID",
    columns = c("GID")
  ) %>% 
    dplyr::distinct(GID) %>% 
    dplyr::mutate(is_SM_gene = TRUE)
  
  
  geneInfo <- dplyr::left_join(x = kmClust, y = geneDesc, by = c("geneId" = "GID")) %>% 
    dplyr::left_join(y = smGenes, by = c("geneId" = "GID"))
  
  
  head(geneInfo)
  
  ## read polII expression data
  expressionData <- get_polII_signal(file = polII_info$polIIExpFile,
                                     title = polII_info$sampleId,
                                     clusterData = geneInfo)
  
  ## add macs2 peak calling information to clusterData
  peakTargets <- peak_target_matrix(sampleInfo = tf_info, position = "best")
  clusterData <- dplyr::left_join(x = expressionData$clusterDf, y = peakTargets, by = "geneId") %>% 
    as.data.frame()
  
  
  matList <- import_profiles(
    exptInfo = sampleInfo,
    geneList = clusterData$geneId,
    source = matrixType,
    # keep = c(0, 0.995),
    up = matrixDim[1], target = matrixDim[2], down = matrixDim[3]
  )
  
  
  ## check the distribution in data
  quantile(matList[[TF_sample]], c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
  
  quantile(matList[[polII_sample]], c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
  
  quantile(matList[[h3_sample]], c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
  
  profileColors[[TF_sample]] <- colorRamp2(
    breaks = quantile(matList[[TF_sample]], c(0.50, 0.995), na.rm = T),
    colors = unlist(strsplit(x = tf_info$color[1], split = ",")))
  
  profileColors[[input_sample]] <- colorRamp2(
    breaks = quantile(matList[[TF_sample]], c(0.50, 0.995), na.rm = T),
    colors = unlist(strsplit(x = input_info$color[1], split = ",")))
  
  profileColors[[polII_sample]] <- colorRamp2(
    breaks = quantile(matList[[polII_sample]], c(0.01, 0.5, 0.995), na.rm = T),
    colors = c("blue", unlist(strsplit(x = polII_info$color[1], split = ","))))
  
  profileColors[[h3_sample]] <- colorRamp2(
    breaks = quantile(matList[[h3_sample]], c(0.2, 0.995), na.rm = T),
    colors = unlist(strsplit(x = hist_info$color[1], split = ",")))
  
  
  profileYlims <- list()
  
  profileYlims[[TF_sample]] <- unname(c(0, quantile(matList[[TF_sample]], 0.996)))
  profileYlims[[input_sample]] <- profileYlims[[TF_sample]]
  
  quantile(expressionData$log2_mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)
  
  polII_color <- colorRamp2(
    breaks = c(0, quantile(expressionData$log2_mat, c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999))),
    colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu"))
  )
  
  peakAnArgs$col[[unname(tfCols$peakType)]] <- structure(
    .Data = RColorBrewer::brewer.pal(
      n = length(na.exclude(unique(clusterData[[tfCols$peakType]]))), name = "Paired"
    ),
    names = na.exclude(unique(clusterData[[tfCols$peakType]]))
  )
  
  ##################################################################################
  ## draw heatmap with all the genes
  cat("## drawing heatmap with all the genes\n")
  
  all_title <- paste0(name, ": all genes", collapse = "")
  
  all_heatmaps <- multi_profile_plots(
    exptInfo = sampleInfo,
    genesToPlot = clusterData$geneId,
    clusters = clusterData,
    profileColors = profileColors,
    matSource = matrixType,
    matBins = matrixDim,
    ylimFraction = profileYlims,
    column_title_gp = gpar(fontsize = 12),
    plotExpression = TRUE,
    expressionData = clusterData,
    expressionColor = polII_color
  )
  
  clusterColors <- all_heatmaps$profileHeatmaps[[TF_sample]]$clusterColor
  
  anGl <- gene_length_heatmap_annotation(
    bedFile = file_genes,
    genes = clusterData$geneId,
    axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"), side = "bottom", labels_rot = 90)
  )
  
  htlist_allGenes <- all_heatmaps$profileHeatmaps[[TF_sample]]$rowGroupHt +
    anGl$an +
    all_heatmaps$profileHeatmaps[[input_sample]]$heatmap +
    all_heatmaps$profileHeatmaps[[TF_sample]]$heatmap +
    all_heatmaps$expressionHeatmaps[[polII_sample]] +
    all_heatmaps$profileHeatmaps[[polII_sample]]$heatmap +
    all_heatmaps$profileHeatmaps[[h3_sample]]$heatmap
  
  # draw Heatmap and add the annotation name decoration
  pdf(file = paste0(outPrefix_all, ".pdf", collapse = ""), width = 16, height = 12)
  htlist_allGenes_draw <- draw(htlist_allGenes,
                               main_heatmap = tf_info$profileName,
                               column_title = all_title,
                               column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                               row_sub_title_side = "left",
                               gap = unit(5, "mm"),
                               padding = unit(rep(0.5, times = 4), "cm")
  )
  
  dev.off()
  
  
  
  ##################################################################################
  ## Draw profile heatmap with the genes for which peak was called by macs2
  
  cat("## drawing heatmap with the genes for which peak was called by macs2\n")
  
  peaks_title <- paste(name, ": macs2 targets", sep = "")
  
  peak_genes <- clusterData[which(clusterData[[unname(tfCols$hasPeak)]]), ]
  rownames(peak_genes) <- peak_genes$geneId
  
  
  peak_heatmaps <- multi_profile_plots(
    exptInfo = sampleInfo,
    genesToPlot = peak_genes$geneId,
    clusters = peak_genes,
    clusterColor = clusterColors,
    profileColors = profileColors,
    matSource = matrixType,
    matBins = matrixDim,
    ylimFraction = profileYlims,
    column_title_gp = gpar(fontsize = 12),
    plotExpression = TRUE,
    expressionData = peak_genes,
    expressionColor = polII_color
  )
  
  
  matData <- as.data.frame(peak_heatmaps$profileHeatmaps$An_kdmB_20h_HA_1$heatmap@matrix) %>% 
    tibble::rownames_to_column(var = "geneId") %>% 
    dplyr::mutate_if(.predicate = is.numeric, .funs = list(~round(., digits = 2))) %>% 
    dplyr::left_join(y = dplyr::select(peak_genes, geneId, cluster, starts_with("peakDist")), by = "geneId") %>% 
    dplyr::select(geneId, cluster, starts_with("peakDist"), everything()) %>% 
    dplyr::arrange(cluster, desc(!!as.name(tfCols$peakDist)))
  
  readr::write_tsv(x = matData, file = paste(outPrefix_peaks, ".matData.tab", sep = ""))
  
  ## peak type annotation
  peakTypeDf <- peak_genes[c(unname(tfCols$peakType))]
  
  peaks_an <- HeatmapAnnotation(
    df = peakTypeDf,
    which = "row",
    col = peakAnArgs$col,
    show_legend = TRUE,
    show_annotation_name = FALSE,
    annotation_width = rep(anWidth, 2),
    gap = unit(2, "mm")
  )
  
  ## peak pval heatmap
  quantile(peak_genes[[unname(tfCols$peakPval)]], c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
  
  peakPvalCol <- colorRamp2(
    breaks = c(1, quantile(peak_genes[[unname(tfCols$peakPval)]], c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 0.95), na.rm = T)),
    colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")))
  
  
  peak_pvalHt <- signal_heatmap(log2_matrix = as.matrix(peak_genes[unname(tfCols$peakPval)]),
                                col_title = "macs2 -log10(pval)",
                                legend_title = "macs2 -log10(pval)",
                                color = peakPvalCol)
  
  
  peaks_gl <- gene_length_heatmap_annotation(
    bedFile = file_genes, genes = peak_genes$geneId,
    axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"), side = "bottom", labels_rot = 90)
  )
  
  
  peaks_htlist <- peak_heatmaps$profileHeatmaps[[TF_sample]]$rowGroupHt +
    peaks_gl$an +
    peak_heatmaps$profileHeatmaps[[input_sample]]$heatmap +
    peak_heatmaps$profileHeatmaps[[TF_sample]]$heatmap +
    peak_pvalHt +
    peaks_an +
    peak_heatmaps$expressionHeatmaps[[polII_sample]] +
    peak_heatmaps$profileHeatmaps[[polII_sample]]$heatmap +
    peak_heatmaps$profileHeatmaps[[h3_sample]]$heatmap
  
  
  
  
  # draw Heatmap and add the annotation name decoration
  pdf(file = paste(outPrefix_peaks, ".pdf", sep = ""), width = 17, height = 12)
  
  draw(peaks_htlist,
       # main_heatmap = tf_info$profileName,
       column_title = peaks_title,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_sub_title_side = "left",
       # heatmap_legend_side = "bottom",
       # annotation_legend_side = "bottom",
       gap = unit(5, "mm"),
       padding = unit(rep(0.5, times = 4), "cm")
  )
  
  
  ## decorate the annotations
  add_annotation_titles(annotations = c(colnames(peakTypeDf)), anTitle = anLables)
  
  dev.off()
  
  
  peaks_htlist2 <- NULL
  peaksSplit <- NULL
  
  ## make sure that the order of genes in the heatmap list and in the dataframe is same
  if(all(rownames(peaks_htlist@ht_list[[tf_info$profileName]]@matrix) == peak_genes$geneId)){
    ## set the row order by decreasing gene length
    # rowOrd_peaks <- order(peak_genes$length, decreasing = TRUE)
    # sliceN <- 1
    
    ## set the row order by peak enrichment
    rowOrd_peaks <- order(peak_genes[[unname(tfCols$peakPval)]], decreasing = TRUE)
    sliceN <- 1
    peaksOrdOutFile <- paste(outPrefix_peaks, "_ordSignal.pdf", sep = "")
    peaksSplit <- rep(1, nrow(peak_genes))
  }
  
  peaks_title <- paste0(name, ": macs2 targets (order: -log10(pVal))", collapse = "")
  
  # draw Heatmap with row order by gene length
  peaks_htlist2 <- peaks_htlist2 +
    peaks_gl$an +
    peak_heatmaps$profileHeatmaps[[input_sample]]$heatmap +
    peak_heatmaps$profileHeatmaps[[TF_sample]]$heatmap +
    peak_pvalHt +
    peaks_an +
    peak_heatmaps$expressionHeatmaps[[polII_sample]] +
    peak_heatmaps$profileHeatmaps[[polII_sample]]$heatmap +
    peak_heatmaps$profileHeatmaps[[h3_sample]]$heatmap
  
  
  
  # draw Heatmap and add the annotation name decoration
  pdf(file = peaksOrdOutFile, width = 17, height = 12)
  
  draw(peaks_htlist2,
       main_heatmap = tf_info$profileName,
       column_title = peaks_title,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_sub_title_side = "left",
       gap = unit(5, "mm"),
       row_order = rowOrd_peaks,
       split = peaksSplit,
       padding = unit(rep(0.5, times = 4), "cm")
  )
  
  
  ## decorate the annotations
  add_annotation_titles(annotations = c(colnames(peakTypeDf)), anTitle = anLables)
  
  dev.off()
  
  
  ## genes ordered by ATG distance from peak
  peaks_htlist3 <- NULL
  
  ## set the row order by distance of the peak from ATG
  if(all(rownames(peaks_htlist@ht_list[[tf_info$profileName]]@matrix) == peak_genes$geneId)){
    rowOrd_peaks <- order(peak_genes[[unname(tfCols$peakDist)]], decreasing = TRUE)
    sliceN <- length(unique(peak_genes$cluster))
    peaks_htlist3 <- peak_heatmaps$profileHeatmaps[[TF_sample]]$rowGroupHt
    peaksOrdOutFile <- paste(outPrefix_peaks, "_ordPeakDist.pdf", sep = "")
  }
  
  peaks_title <- paste0(name, ": macs2 targets (order: peak distance)", collapse = "")
  
  # draw Heatmap with row order by gene length
  peaks_htlist3 <- peaks_htlist3 +
    peaks_gl$an +
    peak_heatmaps$profileHeatmaps[[input_sample]]$heatmap +
    peak_heatmaps$profileHeatmaps[[TF_sample]]$heatmap +
    peak_pvalHt +
    peaks_an +
    peak_heatmaps$expressionHeatmaps[[polII_sample]] +
    peak_heatmaps$profileHeatmaps[[polII_sample]]$heatmap +
    peak_heatmaps$profileHeatmaps[[h3_sample]]$heatmap
  
  
  
  # draw Heatmap and add the annotation name decoration
  pdf(file = peaksOrdOutFile, width = 17, height = 12)
  
  draw(peaks_htlist3,
       main_heatmap = tf_info$profileName,
       column_title = peaks_title,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_sub_title_side = "left",
       gap = unit(5, "mm"),
       row_order = rowOrd_peaks,
       padding = unit(rep(0.5, times = 4), "cm")
  )
  
  
  ## decorate the annotations
  add_annotation_titles(annotations = c(colnames(peakTypeDf)), anTitle = anLables)
  
  dev.off()
  
  ##################################################################################
  
  ## Draw heatmap with only top 10% expressed genes
  
  cat("## drawing heatmap with top 10% expressed genes\n")
  
  expressed_title <- paste0(name, ": top 10% expressed polII genes", collapse = "")
  
  ## extract top 10% expressed genes from polII data
  expressed_genes <- clusterData[which(clusterData[[unname(polIICols$is_expressed)]]), ]
  
  
  expressed_heatmaps <- multi_profile_plots(
    exptInfo = sampleInfo,
    genesToPlot = expressed_genes$geneId,
    clusters = expressed_genes,
    clusterColor = clusterColors,
    profileColors = profileColors,
    matSource = matrixType,
    matBins = matrixDim,
    ylimFraction = profileYlims,
    plotExpression = TRUE,
    expressionData = expressed_genes,
    column_title_gp = gpar(fontsize = 12),
    expressionColor = polII_color
  )
  
  
  exp_gl <- gene_length_heatmap_annotation(
    bedFile = file_genes, genes = expressed_genes$geneId,
    axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"), side = "bottom", labels_rot = 90)
  )
  
  
  expressed_htlist <- expressed_heatmaps$profileHeatmaps[[TF_sample]]$rowGroupHt +
    exp_gl$an +
    expressed_heatmaps$profileHeatmaps[[input_sample]]$heatmap +
    expressed_heatmaps$profileHeatmaps[[TF_sample]]$heatmap +
    expressed_heatmaps$expressionHeatmaps[[polII_sample]] +
    expressed_heatmaps$profileHeatmaps[[polII_sample]]$heatmap +
    expressed_heatmaps$profileHeatmaps[[h3_sample]]$heatmap
  
  
  ## row order as decreasing polII signal
  if(all(rownames(expressed_htlist@ht_list[[tf_info$profileName]]@matrix) == expressed_genes$geneId)){
    expRowOrd <- order(expressed_genes[[polII_sample]], decreasing = TRUE)
  }
  
  
  ## draw Heatmap and add the annotation name decoration
  pdf(file = paste0(outPrefix_expressed, ".pdf", collapse = ""), width = 15, height = 12)
  
  draw(expressed_htlist,
       main_heatmap = tf_info$profileName,
       column_title = expressed_title,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_sub_title_side = "left",
       gap = unit(5, "mm"),
       row_order = expRowOrd,
       padding = unit(rep(0.5, times = 4), "cm")
  )
  
  dev.off()
  
  
  ##################################################################################
  ## Draw profile heatmap with SM genes
  
  cat("## drawing heatmap with SM genes\n")
  
  sm_title <- paste0(name, ": Secondary metabolites genes", collapse = "")
  
  sm_genes <- clusterData[which(clusterData$is_SM_gene), ]
  
  
  sm_heatmaps <- multi_profile_plots(
    exptInfo = sampleInfo,
    genesToPlot = sm_genes$geneId,
    clusters = sm_genes,
    clusterColor = clusterColors,
    profileColors = profileColors,
    matSource = matrixType,
    matBins = matrixDim,
    ylimFraction = profileYlims,
    plotExpression = TRUE,
    expressionData = sm_genes,
    column_title_gp = gpar(fontsize = 12),
    expressionColor = polII_color
  )
  
  
  sm_gl <- gene_length_heatmap_annotation(
    bedFile = file_genes, genes = sm_genes$geneId,
    axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"), side = "bottom", labels_rot = 90)
  )
  
  
  sm_htlist <- sm_heatmaps$profileHeatmaps[[TF_sample]]$rowGroupHt +
    sm_gl$an +
    sm_heatmaps$profileHeatmaps[[input_sample]]$heatmap +
    sm_heatmaps$profileHeatmaps[[TF_sample]]$heatmap +
    sm_heatmaps$expressionHeatmaps[[polII_sample]] +
    sm_heatmaps$profileHeatmaps[[polII_sample]]$heatmap +
    sm_heatmaps$profileHeatmaps[[h3_sample]]$heatmap
  
  
  # draw Heatmap and add the annotation name decoration
  pdf(file = paste0(outPrefix_sm, ".pdf", collapse = ""), width = 17, height = 12)
  
  draw(sm_htlist,
       main_heatmap = tf_info$profileName,
       # annotation_legend_list = list(sm_profile$legend),
       column_title = sm_title,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_sub_title_side = "left",
       gap = unit(5, "mm"),
       padding = unit(rep(0.5, times = 4), "cm")
  )
  
  dev.off()
  
  ##################################################################################
  
  ## plot polII expressed and TF bound genes together
  
  cat("## drawing heatmap with top 10% polII expressed and TF bound genes together\n")
  
  pkExp_title <- paste0(name, ": polII expressed and TF bound genes", collapse = "")
  
  
  ## select the genes which are expressed in polII sample OR have TSS peak
  pkExp_genes <- filter_at(.tbl = clusterData, .vars = unname(c(polIICols$is_expressed, tfCols$hasPeak)),
                           .vars_predicate = any_vars(. == "TRUE"))
  
  
  pkExp_heatmaps <- multi_profile_plots(
    exptInfo = sampleInfo,
    genesToPlot = pkExp_genes$geneId,
    clusters = pkExp_genes,
    clusterColor = clusterColors,
    profileColors = profileColors,
    matSource = matrixType,
    matBins = matrixDim,
    ylimFraction = profileYlims,
    plotExpression = TRUE,
    expressionData = pkExp_genes,
    column_title_gp = gpar(fontsize = 12),
    expressionColor = polII_color
  )
  
  
  pkExp_gl <- gene_length_heatmap_annotation(
    bedFile = file_genes, genes = pkExp_genes$geneId,
    axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"), side = "bottom", labels_rot = 90)
  )
  
  ## annotations
  peakAnArgs$col[["is_SM_gene"]] <- c("TRUE" = "orange", "FALSE" = "grey95")
  peakAnArgs$col[["is_TF"]] <- c("TRUE" = "blue", "FALSE" = "grey95")
  peakAnArgs$col[[unname(polIICols$is_expressed)]] <- structure(c("green4", "grey95"), names = c("TRUE", "FALSE"))
  
  pkExp_An <- HeatmapAnnotation(
    df = pkExp_genes[c("is_SM_gene", unname(polIICols$is_expressed), unname(tfCols$peakType), unname(tfCols$tesPeakType))],
    which = "row",
    col = peakAnArgs$col,
    show_annotation_name = FALSE,
    show_legend = TRUE,
    annotation_width = rep(anWidth, 5),
    gap = unit(2, "mm")
  )
  
  
  ## heatmap list
  pkExp_htlist <- pkExp_heatmaps$profileHeatmaps[[TF_sample]]$rowGroupHt +
    pkExp_gl$an +
    pkExp_heatmaps$profileHeatmaps[[input_sample]]$heatmap +
    pkExp_heatmaps$profileHeatmaps[[TF_sample]]$heatmap +
    pkExp_heatmaps$expressionHeatmaps[[polII_sample]] +
    pkExp_heatmaps$profileHeatmaps[[polII_sample]]$heatmap +
    pkExp_An +
    pkExp_heatmaps$profileHeatmaps[[h3_sample]]$heatmap
  
  
  ## row order as decreasing polII signal
  if(all(rownames(pkExp_htlist@ht_list[[tf_info$profileName]]@matrix) == pkExp_genes$geneId)){
    pkExpOrd <- order(pkExp_genes[[unname(tfCols$peakDist)]], pkExp_genes[[polII_sample]], decreasing = TRUE)
  }
  
  
  # draw Heatmap and add the annotation name decoration
  pdf(file = paste0(outPrefix_pkExp, ".pdf", collapse = ""), width = 17, height = 12)
  
  draw(pkExp_htlist,
       main_heatmap = tf_info$profileName,
       column_title = pkExp_title,
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       row_sub_title_side = "left",
       row_order = pkExpOrd,
       gap = unit(5, "mm"),
       padding = unit(rep(0.5, times = 4), "cm")
  )
  
  
  ## decorate the annotations
  add_annotation_titles(
    annotations = c("is_SM_gene", unname(polIICols$is_expressed), unname(tfCols$peakType)),
    anTitle = anLables)
  
  dev.off()
  
  
  ##################################################################################
  ## store the clusterData
  write.table(clusterData, file = tf_info$mergedDataFile,
              col.names = T, row.names = F, sep = "\t", quote = F)
  
  ## write data to excel file
  wb <- openxlsx::createWorkbook(creator = "Lakhansing Pardehi")
  openxlsx::addWorksheet(wb = wb, sheetName = name)
  openxlsx::writeData(
    wb = wb, sheet = 1, x = clusterData,
    startCol = 1, startRow = 1, withFilter = TRUE,
    keepNA = TRUE, na.string = "NA"
  )
  headerStyle <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#e6e6e6")
  openxlsx::addStyle(wb = wb, sheet = 1, style = headerStyle, rows = 1, cols = 1:ncol(clusterData))
  openxlsx::setColWidths(wb = wb, sheet = 1, cols = 1, widths = "auto")
  openxlsx::freezePane(wb = wb, sheet = 1, firstActiveRow = 2, firstActiveCol = 2)
  
  # openxlsx::openXL(wb)
  openxlsx::saveWorkbook(wb = wb, file = excelOut, overwrite = TRUE)
  
  ##################################################################################
  
  
  cat(TF_sample, ": Done\n")
  
  
}


stopCluster(cl)




## all the functions required for integrating multiple TF data


##################################################################################
## generate experiment data
getSampleInformation = function(exptInfoFile, samples, polII_path, TF_path){
  
  ## read the experiment sample details and select only those which are to be plotted
  exptData = fread(input = exptInfoFile, 
                   sep = "\t",
                   stringsAsFactors = F, 
                   header = T,
                   data.table = F) %>% 
    filter(Sample_ID %in% samples)
  
  exptData$Sample_ID = factor(exptData$Sample_ID, levels = samples)
  
  exptData = exptData[order(exptData$Sample_ID), ] %>%
    mutate_if(is.factor, as.character) %>%
    mutate(
      profileName = paste(Sample_ID, "profile", sep = "_"),
      matFile = if_else(
        IP_tag == "polII", 
        paste(polII_path, "/", Sample_ID, "/", Sample_ID, "_normalized_profile.tab.gz", sep = ""),
        paste(TF_path, "/", Sample_ID, "/", Sample_ID, "_normalized_profile.tab.gz", sep = "")
      ),
      polIIExpFile = if_else(
        IP_tag == "polII",
        paste(polII_path, "/", Sample_ID, "/", Sample_ID, "_normalizedExpression.tab", sep = ""),
        "NA"
      ),
      polIIExpMat = if_else(
        IP_tag == "polII",
        paste(polII_path, "/", Sample_ID, "/", Sample_ID, "_polii_expr.tab.rel.mat", sep = ""),
        "NA"
      ),
      clusterFile = if_else(
        IP_tag == "polII",
        "NA",
        paste(TF_path, "/", Sample_ID, "/", Sample_ID, ".kmeans.clusters.txt", sep = "")
      ),
      tfPeakFile = if_else(
        IP_tag == "polII",
        "NA",
        paste(TF_path, "/", Sample_ID, "/", Sample_ID, "_allGenes_clusters.tab", sep = "")
      ),
      narrowpeakFile = if_else(
        IP_tag == "polII",
        "NA",
        paste(TF_path, "/", Sample_ID, "/", Sample_ID, "_narrow_peaks.narrowPeak", sep = "")
      ),
      narrowpeakAnno = if_else(
        IP_tag == "polII",
        "NA",
        paste(TF_path, "/", Sample_ID, "/", Sample_ID, "_narrow_summits_nearest_ATG.tab", sep = "")
      ),
      broadpeakFile = if_else(
        IP_tag == "polII",
        "NA",
        paste(TF_path, "/", Sample_ID, "/", Sample_ID, "_broad_peaks.broadPeak", sep = "")
      ),
      broadpeakAnno = if_else(
        IP_tag == "polII",
        "NA",
        paste(TF_path, "/", Sample_ID, "/", Sample_ID, "_broad_center_nearest_ATG.tab", sep = "")
      )     
    ) %>%
    mutate_at(c("polIIExpFile", "polIIExpMat", "clusterFile", "tfPeakFile", "narrowpeakFile", "narrowpeakAnno", "broadpeakFile", "broadpeakAnno"), 
              .funs = funs(na_if(., "NA")))
  
  
  return(exptData)
}

##################################################################################




##################################################################################
## function to plot the gene information as Heatmap Annotation
getGeneInfo = function(dataFile, orderFile){
  ## arguments
  ## dataFile:  gene info data file
  ## orderFile: SM cluster order file
  
  ## get geneInfo data
  geneData = fread(input = dataFile, header = T, drop = c(2,3),  stringsAsFactors = F, sep = "\t", data.table = F) %>%
    # filter(is_SM_gene == "TRUE") %>%
    select(gene, is_SM_gene, SM_cluster, is_TF)
  
  
  return(geneData)
}


## get the polII expression values and top 10% status
get_polII_expressions = function(genesDf, exptInfo){
  
  if(is.null(exptInfo$polIIExpFile)){
    warning("Expression data not found...")
    return(genesDf)
  }
  
  for(i in 1:nrow(exptInfo)){
    
    if(exptInfo$IP_tag[i] == "polII" & !is.na(exptInfo$polIIExpFile[i])){
      
      df = fread(input = exptInfo$polIIExpFile[i], header = T, stringsAsFactors = F, sep = "\t", data.table = F)
      genesDf = right_join(x = genesDf, y = df, by = c("gene" = "gene"))
      
    }
  }
  
  return(genesDf)
}

##################################################################################




##################################################################################
## get the TF peak data for each gene
get_TF_binding_data = function(TF_path, genesDf, exptInfo, peakCol){
  
  if (is.null(exptInfo$tfPeakFile)) {
    warning("TF binding data not found...")
    return(genesDf)
  }
  
  for(i in 1:nrow(exptInfo)){
    if(exptInfo$IP_tag[i] != "polII" & !is.na(exptInfo$tfPeakFile[i])){
      
      colName = paste0(peakCol, "(", exptInfo$Sample_ID[i], ")", collapse = "")
      
      df = fread(input = exptInfo$tfPeakFile[i], header = T, stringsAsFactors = F, sep = "\t", data.table = F) %>%
        select(gene, !!colName := !! peakCol)
      
      genesDf = left_join(x = genesDf, y = df, by = c("gene" = "gene"))
      
    }
  }
  
  return(genesDf)
}

##################################################################################



##################################################################################

## generate the multi sample comparison plot for SM genes
SM_genesComparisonPlot = function(data, orderFile, exptInfo, outPath, name){
  
  plotTitle = paste("Secondary metabolite genes for comparison", name, "\n")
  outFile = paste0(outPath, "/", name, "_compare_SM_genes.png", collapse = "")
  
  ## prepare SM dataframe
  smData = data %>% filter(is_SM_gene == TRUE)
  
  clusterOrder = fread(input = orderFile, header = T, stringsAsFactors = F, data.table = F)
  smData$SM_cluster = factor(smData$SM_cluster, levels = clusterOrder$clusters)
  smData = smData[order(smData$SM_cluster), ]
  
  row.names(smData) = smData$Gene
  
  
  htSplit = smData["SM_cluster"]
  
  
  ## clusterAnnotation
  clusterIds = as.character(unique(smData$SM_cluster))
  cluster_colors = structure(rep(x = c("lightgrey", "darkgrey"), length.out = length(clusterIds)), names = clusterIds)
  
  clustAn = rowAnnotation(df = smData["SM_cluster"],
                          col = list(SM_cluster = cluster_colors),
                          width = unit(1.5, "cm"),
                          show_legend = FALSE
  )
  
  
  ## TF annotation
  tfAn = rowAnnotation(df = smData["is_TF"],
                       col = list(is_TF = c("TRUE" = "blue", "FALSE" = "grey95")),
                       width = unit(1.5, "cm"),
                       show_legend = FALSE
  )
  
  
  ## plot polII data
  polII_ids = exptInfo$Sample_ID[which(exptInfo$IP_tag == "polII")]
  
  polII_mat = as.matrix(smData[polII_ids])
  
  polII_log2_mat = log2(polII_mat + 1)
  
  quantile(polII_log2_mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)
  
  
  polII_color = colorRamp2(breaks = seq(quantile(polII_log2_mat, 0.5), quantile(polII_log2_mat, 0.995), length.out = 9), colors = c(colors_seq9$Reds))
  
  
  polII_ht = Heatmap(matrix = polII_log2_mat,
                     name = "log2(polII_expression)",
                     col = polII_color,
                     split = htSplit,
                     cluster_rows = TRUE, 
                     cluster_columns = FALSE,
                     show_row_names = TRUE,
                     column_names_gp = gpar(fontsize = 16),
                     row_names_side = "right",
                     row_names_gp = gpar(fontsize = 3),
                     row_title_rot = 0,
                     row_title_gp = gpar(fontsize = 12),
                     show_heatmap_legend = FALSE,
                     heatmap_legend_param = list(title = "log2(polII_expression)\n",
                                                 title_gp = gpar(fontsize = 18, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 14),
                                                 legend_height = unit(4, "cm")
                                                 # legend_width = unit(2, "cm")
                     )
  )
  
  
  htLgd = Legend(title = "log2(polII_expression)",
                 at = as.numeric(sprintf("%.0f", attributes(polII_color)$breaks)),
                 # labels = sprintf("%.0f", attributes(polII_color)$breaks),
                 col_fun = polII_color,
                 type = "grid",
                 direction = "horizontal",
                 title_gp = gpar(fontsize = 12, fontface = "bold"),
                 labels_gp = gpar(fontsize = 8),
                 legend_width = unit(4, "cm"),
                 grid_height = unit(5, "mm")
                 # grid_width = unit(5, "mm")
  )
  
  
  ## is_expressed(polII_sample) annotations
  isExpressedDf = smData[paste("is_expressed(", polII_ids, ")", sep = "")]
  
  expressColor = sapply(X = names(isExpressedDf), 
                        FUN = function(x){structure(c("green4", "grey95"), names = c("TRUE", "FALSE"))}, 
                        simplify = F, 
                        USE.NAMES = T
  )
  
  expressedAnn = HeatmapAnnotation(df = isExpressedDf, 
                                   which = "row", 
                                   col = expressColor,
                                   show_legend = F, 
                                   annotation_width = unit(rep(1.5, times = length(names(isExpressedDf))), "cm"),
                                   gap = unit(3, "mm")
  )
  
  
  ## has_TSS_peak(TF_sample) annotations
  tfSampleIds = exptInfo$Sample_ID[which(exptInfo$IP_tag != "polII")]
  
  hasPeakDf = smData[paste("has_TSS_peak(", tfSampleIds, ")", sep = "")]
  
  tfColor = sapply(X = names(hasPeakDf),
                   FUN = function(x){structure(c("sienna", "grey95"), names = c("TRUE", "FALSE"))},
                   simplify = F,
                   USE.NAMES = T
  )
  
  hasPeakAnn = HeatmapAnnotation(df = hasPeakDf, 
                                 which = "row", 
                                 col = tfColor, 
                                 show_legend = F,
                                 annotation_width = unit(rep(1.5, times = length(names(hasPeakDf))), "cm"),
                                 gap = unit(3, "mm")
  )
  
  
  annLgd = Legend(title = "\nGene categories",
                  at = c("Transcription Factor", "Expressed gene", "TF binding peak at TSS"),
                  type = "points",
                  pch = 15, 
                  size = unit(5, "mm"),
                  grid_height = unit(5, "mm"),
                  grid_width = unit(5, "mm"),
                  legend_gp = gpar(col = c("blue", "green4", "sienna")),
                  title_gp = gpar(fontsize = 12, fontface = "bold"),
                  labels_gp = gpar(fontsize = 8)
  )
  
  
  
  ## heatmap drawing
  htList = clustAn + tfAn + polII_ht + expressedAnn + hasPeakAnn
  
  imageWidth = (length(polII_ids) * 4000) + (length(tfSampleIds) * 2000) + 4000
  png(filename = outFile, width = imageWidth, height = 30000, res = 1200)
  
  ## draw Heatmap list
  draw(htList,
       column_title = plotTitle,
       column_title_gp = gpar(fontsize = 22, fontface = "bold"),
       padding = unit(c(3, 1, 1, 1), "cm"),
       row_sub_title_side = "left",
       annotation_legend_list = list(htLgd, annLgd)
  )
  
  
  ## SM_cluster annotation name
  decorate_annotation(annotation = "SM_cluster", 
                      code = {
                        grid.text(label = "SM clusters", x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm") ,
                                  default.units = "npc", just = "right", rot = 90,
                                  gp = gpar(fontsize = 16))
                      }, 
                      slice = length(clusterIds)
  )
  
  
  ## is_TF annotation name
  decorate_annotation(annotation = "is_TF", 
                      code = {
                        grid.text(label = "Transcription Factor", 
                                  x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm") , 
                                  default.units = "npc", just = "right", rot = 90,
                                  gp = gpar(fontsize = 16))
                      },
                      slice = length(clusterIds)
  )
  
  
  ## isExpressed annotation names
  for(an in colnames(isExpressedDf)){
    txt = gsub("is_expressed", "is expressed\n", an) %>% gsub("\\(|\\)", "", .)
    
    decorate_annotation(annotation = an, 
                        code = {
                          grid.text(label = txt, x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                    default.units = "npc", just = "right", rot = 90,
                                    gp = gpar(fontsize = 16))
                        },
                        slice = length(clusterIds)
    )
  }
  
  
  ## hasPeak annotation names
  for(an in colnames(hasPeakDf)){
    txt = gsub("has_TSS_peak", "has TSS peak\n", an) %>% gsub("\\(|\\)", "", .)
    
    decorate_annotation(annotation = an, 
                        code = {
                          grid.text(label = txt, x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                    default.units = "npc", just = "right", rot = 90,
                                    gp = gpar(fontsize = 16))
                        },
                        slice = length(clusterIds)
    )
  }
  
  
  dev.off()
  
  
}


##################################################################################



##################################################################################

## for a set of genes, generate a multi TF heatmap
geneset_multiTF_plot = function(data, exptInfo, outPath, outPrefix, plotTitle, name, res, numClust){
  
  outFile = paste0(outPath, "/", outPrefix, ".png", collapse = "")
  
  ## empty heatmap list
  htList = NULL
  lgdList = list()
  
  
  polII_ids = exptInfo$Sample_ID[which(exptInfo$IP_tag == "polII")]
  polII_expIds = paste("is_expressed(", polII_ids, ")", sep = "")
  
  ## is_expressed(polII_sample) annotations
  expressColor = sapply(X = polII_expIds, 
                        FUN = function(x){structure(c("green4", "grey95"), names = c("TRUE", "FALSE"))}, 
                        simplify = F, 
                        USE.NAMES = T
  )
  
  expressedAn = HeatmapAnnotation(df = data[polII_expIds], 
                                  which = "row", 
                                  col = expressColor,
                                  show_legend = FALSE,
                                  annotation_width = unit(rep(1, times = length(polII_expIds)), "cm"),
                                  gap = unit(3, "mm")
  )
  
  
  
  ## SM gene annotation
  smGeneAn = HeatmapAnnotation(df = data["is_SM_gene"],
                               which = "row",
                               col = list(is_SM_gene = c("TRUE" = "orange", "FALSE" = "grey95")),
                               annotation_width = unit(1, "cm"),
                               show_legend = FALSE
  )
  
  
  ## TF annotation
  tfAn = rowAnnotation(df = data["is_TF"],
                       col = list(is_TF = c("TRUE" = "blue", "FALSE" = "grey95")),
                       annotation_width = unit(1, "cm"),
                       show_legend = FALSE
  )
  
  
  
  ## clustering polII expression data
  polII_mat = as.matrix(data[polII_ids])
  
  polII_log2_mat = log2(polII_mat + 1)
  dend = hclust(dist(polII_log2_mat), method = "ward.D")
  
  
  ## polII Heatmap
  quantile(polII_log2_mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)
  
  
  polII_color = colorRamp2(breaks = seq(quantile(polII_log2_mat, 0), quantile(polII_log2_mat, 0.99), length.out = 9), colors = c(colors_seq9$Reds))
  # colorRamp2(breaks = quantile(polII_log2_mat, 0, 0.99),colors =  c("white", "red"))
  
  polII_ht = Heatmap(matrix = polII_log2_mat,
                     name = name,
                     col = polII_color,
                     cluster_rows = dend,
                     row_dend_width = unit(3, "cm"),
                     row_dend_reorder = TRUE,
                     row_dend_side = "left",
                     cluster_columns = FALSE,
                     show_row_names = FALSE,
                     column_names_gp = gpar(fontsize = 10),
                     row_title_rot = 0,
                     row_title_gp = gpar(fontsize = 12),
                     show_heatmap_legend = FALSE,
                     heatmap_legend_param = list(title = "log2(polII_expression)",
                                                 title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 8),
                                                 color_bar = "continuous", 
                                                 legend_direction = "horizontal", 
                                                 legend_width = unit(4, "cm"), 
                                                 # legend_height = unit(4, "cm"),
                                                 grid_height = unit(5, "mm")
                                                 # grid_width = unit(5, "mm")
                     ),
                     width = 1 * ncol(polII_log2_mat)
  )
  
  ## heatmap legend
  htLgd = Legend(title = "log2(polII_expression)",
                 at = as.numeric(sprintf("%.0f", attributes(polII_color)$breaks)),
                 # labels = sprintf("%.0f", attributes(polII_color)$breaks),
                 col_fun = polII_color,
                 type = "grid",
                 direction = "horizontal",
                 title_gp = gpar(fontsize = 12, fontface = "bold"),
                 labels_gp = gpar(fontsize = 8),
                 legend_width = unit(4, "cm"),
                 grid_height = unit(5, "mm")
                 # grid_width = unit(5, "mm")
  )
  
  
  ## has_TSS_peak(TF_sample) annotations
  tfSampleIds = exptInfo$Sample_ID[which(exptInfo$IP_tag != "polII")]
  
  hasPeakDf = data[paste("has_TSS_peak(", tfSampleIds, ")", sep = "")]
  
  tfColor = sapply(X = names(hasPeakDf),
                   FUN = function(x){structure(c("sienna", "grey95"), names = c("TRUE", "FALSE"))},
                   simplify = F,
                   USE.NAMES = T
  )
  
  hasPeakAn = HeatmapAnnotation(df = hasPeakDf, 
                                which = "row", 
                                col = tfColor, 
                                show_legend = F,
                                annotation_width = unit(rep(1, times = length(names(hasPeakDf))), "cm"),
                                gap = unit(3, "mm")
  )
  
  
  ## legend for all the annotations
  annLgd = Legend(title = "\nGene categories",
                  at = c("SM gene", "Transcription Factor", "Expressed gene", "TF binding peak at TSS"),
                  type = "points",
                  pch = 15, 
                  size = unit(5, "mm"),
                  grid_height = unit(5, "mm"),
                  grid_width = unit(5, "mm"),
                  legend_gp = gpar(col = c("orange", "blue", "green4", "sienna")),
                  title_gp = gpar(fontsize = 12, fontface = "bold"),
                  labels_gp = gpar(fontsize = 8)
  )
  
  
  ## cut dendrogram and draw cluster annotation only if number of cuts given
  clustLgd = NULL
  clusterData = NULL
  
  if(!missing(numClust)){
    
    clusterCut = cutree(dend, numClust)
    clusterData = data.frame(gene = names(clusterCut), 
                             cluster = paste("cluster", clusterCut, sep = "_"), 
                             stringsAsFactors = F)
    
    clusterColor <- structure(colors_qual$set1_9[1:numClust], names = unique(clusterData$cluster))
    
    # cluster annotation
    clustAn <- HeatmapAnnotation(df = clusterData["cluster"], 
                                 which = "row",
                                 width = unit(2, "mm"),
                                 col = list(cluster = clusterColor),
                                 show_legend = FALSE
    )
    
    htList = htList + clustAn
    
    clustLgd = Legend(title = "Gene clusters",
                      at = names(clusterColor),
                      type = "points",
                      pch = 15, 
                      size = unit(5, "mm"),
                      grid_height = unit(5, "mm"),
                      grid_width = unit(5, "mm"),
                      legend_gp = gpar(col = clusterColor),
                      title_gp = gpar(fontsize = 12, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8),
                      ncol = 2
    )
    
    data = left_join(x = data, y = clusterData, by = c("gene" = "gene"))
  }
  
  ## all plotting Heatmap and HeatmapAnnotation objects
  htList = htList + smGeneAn + tfAn + polII_ht + expressedAn + hasPeakAn
  
  lgdList = list(htLgd, annLgd)
  
  if(!missing(numClust)){
    lgdList = list(clustLgd, htLgd, annLgd)
  }
  
  
  ## generate plot
  
  imageWidth = (length(polII_ids) * 750) + (length(tfSampleIds) * 500) + 1500
  imageHeight = nrow(polII_log2_mat) * 5
  imageRes = imageWidth/10
  
  if(!missing(res)){
    imageRes = res
  }
  
  
  png(filename = outFile, width = imageWidth, height = imageHeight, res = imageRes)
  
  ## draw Heatmap list
  draw(htList,
       column_title = plotTitle,
       column_title_gp = gpar(fontsize = 16, fontface = "bold"),
       padding = unit(c(2, 1, 1, 1), "cm"),
       row_dend_side = "left",
       row_sub_title_side = "left",
       annotation_legend_list = lgdList
  )
  
  
  ## SM_cluster annotation name
  decorate_annotation(annotation = "is_SM_gene", 
                      code = {
                        grid.text(label = "SM gene", x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm") ,
                                  default.units = "npc", just = "right", rot = 90,
                                  gp = gpar(fontsize = 10)
                        )
                      }
  )
  
  
  ## is_TF annotation name
  decorate_annotation(annotation = "is_TF", 
                      code = {
                        grid.text(label = "Transcription Factor", 
                                  x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm") , 
                                  default.units = "npc", just = "right", rot = 90,
                                  gp = gpar(fontsize = 10)
                        )
                      }
  )
  
  
  ## isExpressed annotation names
  for(an in polII_expIds){
    txt = gsub("is_expressed", "is expressed\n", an) %>% gsub("\\(|\\)", "", .)
    
    decorate_annotation(annotation = an, 
                        code = {
                          grid.text(label = txt, x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                    default.units = "npc", just = "right", rot = 90,
                                    gp = gpar(fontsize = 10)
                          )
                        }
    )
  }
  
  
  ## hasPeak annotation names
  for(an in colnames(hasPeakDf)){
    txt = gsub("has_TSS_peak", "has TSS peak\n", an) %>% gsub("\\(|\\)", "", .)
    
    decorate_annotation(annotation = an, 
                        code = {
                          grid.text(label = txt, x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                    default.units = "npc", just = "right", rot = 90,
                                    gp = gpar(fontsize = 10)
                          )
                        }
    )
  }
  
  
  dev.off()
  
  return(list(
    "df" = data,
    "dend" = dend
  ))
  
}

##################################################################################




##################################################################################
## function to get the profile plots in list() using EnrichedHeatmap package
getProfilePlotList = function(df, cluster, clusterColor, geneList, colorList = NULL, matrixSource = "deeptools", matrixBins){
  allPlots = list()
  
  ## generate plot for each sample
  for(i in 1:nrow(df)){
    sampleName = df$Sample_ID[i]
    
    cat("Generating profile heatmap for:", sampleName, "...\n")
    
    mat = matrix()
    
    ## read the profile matrix
    if(matrixSource == "deeptools"){
      mat = getDeeptoolsProfileMatrix(file = df$matFile[i], 
                                      signalName = sampleName, 
                                      selectGenes = geneList)
    }
    else if(matrixSource == "miao"){

      mat = getMiaoProfileMatrix(file = df$matFile[i], 
                                 signalName = sampleName, 
                                 selectGenes = geneList,
                                 up = matrixBins[1],
                                 target = matrixBins[2],
                                 down = matrixBins[3])
    }
    
    
    
    
    ## color: generate if colorList = NULL
    tfCol = colorRamp2(quantile(mat, c(0.50, 0.99), na.rm = T), c("white", "red"))
    polIICol = colorRamp2(quantile(mat, c(0.01, 0.5, 0.99), na.rm = T), c("blue", "white", "red"))
    
    colFun = tfColFun
    if(df$IP_tag[i] == "polII"){
      colFun = polIIColFun
    }
    
    ## if the colors are provided with colorList argument, use those colors instead
    if(!is.null(colorList)){
      colFun = colorList[[sampleName]]
    }
    
    ## build the profile heatmap
    profile2 = secondaryProfileHeatmap(profileMat = mat, 
                                       signalName = df$profileName[i], 
                                       clusterData = cluster,
                                       heatmapColor = colFun,
                                       clusterColor = clusterColor)
    
    ## append to the plot list
    allPlots[[ sampleName ]] = profile2
    
    cat("Done!!!\n")
  }
  
  return(allPlots)
}

##################################################################################




##################################################################################

## get the pol-II expression heatmap in list() of class ComplexHeatmap for many samples 
getExpressionHeatmapList = function(expDf, exptInfo, htColor = NULL){
  
  polII_ids = exptInfo$Sample_ID[which(exptInfo$IP_tag == "polII")]
  polII_expIds = paste("is_expressed(", polII_ids, ")", sep = "")
  
  
  ## clustering polII expression data
  rownames(expDf) = expDf$gene
  polII_mat = as.matrix(expDf[polII_ids])
  
  polII_log2_mat = log2(polII_mat + 1)
  
  ## polII Heatmap
  quantile(polII_log2_mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 1), na.rm = T)
  
  
  polII_color = colorRamp2(breaks = seq(quantile(polII_log2_mat, 0), quantile(polII_log2_mat, 0.99), length.out = 9), colors = c(colors_seq9$Reds))
  # colorRamp2(breaks = quantile(polII_log2_mat, 0, 0.99),colors =  c("white", "red"))
  
  if(!is.null(htColor)){
    polII_color = htColor
  }
  
  allPlots = list()
  firstPlot = TRUE
  
  ## generate heatmaps
  for(id in polII_ids){
    ht = polII_heatmap_secondary(log2_matrix = polII_log2_mat[, id], 
                                 col_title = id,
                                 legend_title = "log2(polII_FPKM + 1)", 
                                 color = polII_color,
                                 showLegend = firstPlot)
    
    firstPlot = FALSE
    
    allPlots[[id]] = ht 
  }
  
  allPlots[["polII_color"]] = polII_color
  
  return(allPlots)
}


##################################################################################



##################################################################################

## function to get a heatmap list of provided TF, polII samples to be used directly for plotting
multiProfilePlots = function(exptInfo, expressionData, genesToPlot, clusters, clusterColor, profileColors = NULL, expressionColor = NULL, plotExpression = FALSE, matSource, matBins){
  
  ## exptInfo = experiment info as data frame with information like sampleID, type, path etc
  ## expressionData = an merged dataframe which has info: clustering, polII expression, TF binding status for each gene
  ## genesToPlot = a data frame with gene of interest
  ## clusters = cluster information in forma of dataframe
  ## matSource = miao / deeptools
  ## matBins = c(up, body, down) bin count
  ## clusterColor = cluster color object
  ## profileColors = a named list of color objects for each of the sample
  ## expressionColor = color object for polII expression heatmap
  
  cat("Generating profile heatmaps...\n")
  
  ## profile heatmap list for all samples except first
  profileList = getProfilePlotList(df = exptInfo,
                                   cluster = clusters,
                                   clusterColor = clusterColor,
                                   geneList = genesToPlot, 
                                   colorList = profileColors,
                                   matrixSource = matSource,
                                   matrixBins = matBins)
  
  cat("Done!!!\n")
  
  expHtList = list()
  expHtList[["polII_color"]] = expressionColor
  
  if(isTRUE(plotExpression)){
    
    cat("Generating expression heatmaps...\n")
    
    ## expression heatmap list
    expHtList = getExpressionHeatmapList(expDf = expressionData, 
                                         exptInfo = exptInfo, 
                                         htColor = expressionColor)
    
    cat("Done!!!\n")
  }
  
  
  cat("Generating heatmap list object...\n")
  
  ## build the heatmap list for plotting
  htList = profileList[[exptInfo$Sample_ID[1]]][["rowAnno"]]
  plotGaps = c(2)
  profileColorList = list()
  
  for(i in 1:nrow(exptInfo)){
    sampleID = exptInfo$Sample_ID[i]
    
    ## add plot elements to the heatmap list
    if(exptInfo$IP_tag[i] == "polII"){
      
      ## polII expression + profile
      if(isTRUE(plotExpression)){
        htList = htList + expHtList[[sampleID]] + profileList[[sampleID]][["heatmap"]]
        plotGaps = append(plotGaps, c(2, 10))
      }
      ## only polII profile
      else{
        htList = htList + profileList[[sampleID]][["heatmap"]]
        plotGaps = append(plotGaps, 10)
      }
      
    }
    ## TF profile
    else{
      htList = htList + profileList[[sampleID]][["heatmap"]]
      plotGaps = append(plotGaps, 10)
    }
    
    ## add profile colors to profileColorList
    profileColorList[[sampleID]] = profileList[[sampleID]][["profileColor"]]
    
  }
  
  plotGaps = head(plotGaps, -1)
  
  cat("Done!!!\n")
  
  return(list(
    "heatmapList" = htList,
    "plotGaps" = plotGaps,
    "profileColors" = profileColorList,
    "expressionColor" = expHtList$polII_color
  ))
  
}

##################################################################################



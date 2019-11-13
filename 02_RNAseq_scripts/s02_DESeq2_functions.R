




###########################################################################
## function to plot volcano plot
volcano_plot <- function(df, title = "volcano plot",
                         fdr_col, lfc_col, fdr_cut = 0.05, lfc_cut = 1,
                         ylimit = 150, xlimit = c(-4, 4), geneOfInterest, showNames = TRUE){
  
  if (is.null(df[["geneName"]])) {
    df$geneName = df$geneId
  }
  
  up_cut <- lfc_cut
  down_cut <- -1 * lfc_cut
  
  df <- df %>% 
    dplyr::mutate(
      category = dplyr::case_when(
        !! sym(fdr_col) < !! fdr_cut & !! sym(lfc_col) >= !! up_cut ~ "Significant Up",
        !! sym(fdr_col) < !! fdr_cut & !! sym(lfc_col) <= !! down_cut ~ "Significant Down",
        !! sym(fdr_col) < !! fdr_cut ~ "Significant",
        TRUE ~ "Non-significant"
      )
    )
  
  
  ## squish the value to limits
  df$log10FDR <- -log10(df[[fdr_col]])
  df$log10FDR <- scales::squish(x = df$log10FDR, range = c(0, ylimit))
  df[[lfc_col]] <- scales::squish(x = df[[lfc_col]], range = xlimit)
  
  
  #Drwa Volcano plot
  p <- ggplot() +
    geom_hline(yintercept = -log10(fdr_cut), color = "black", linetype = "dashed") +
    geom_vline(xintercept = -lfc_cut, color = "black", linetype = "dashed") +
    geom_vline(xintercept = lfc_cut, color = "black", linetype = "dashed") +
    # geom_point(mapping = aes(color = significant), alpha=0.6, size=1.75) +
    geom_point(
      data = df,
      mapping = aes(x = !! sym(lfc_col), y = log10FDR, color = category), alpha=0.6, size=1.75
    ) +
    scale_color_manual(
      values = c("Significant Up" = "red", "Significant Down" = "green",
                 "Non-significant" = "black", "Significant" = "grey"), 
      name = "Significance") +
    scale_x_continuous(name = "log2(fold change)", limits = xlimit, expand = expand_scale(mult = 0.02)) +
    scale_y_continuous(name = "-log10(q-value)", limits = c(0, ylimit), expand = expand_scale(mult = 0.02)) +
    guides(color = guide_legend(nrow = 2, byrow = T)) +
    theme_bw() +
    theme(legend.background = element_rect(colour = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(face = "bold", size = 16),
          legend.position = "bottom",
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.margin = unit(c(1,1,1,1),"cm"),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 16)) +
    ggtitle(title)
  
  
  ## highlight genes of interest
  if(!missing(geneOfInterest)){
    tmpDf <- dplyr::filter(df, geneName %in% geneOfInterest)
    
    ## draw the points
    p <- p +
      geom_point(data = tmpDf,
                 color = "black",
                 shape = 1
      )
    
    ## show the gene lables
    if(isTRUE(showNames)){
      p <- p +
        geom_text_repel(data = tmpDf, mapping = aes(label = geneName),
                        segment.color = '#cccccc',
                        segment.size = 1,
                        size = 5) 
    }
  }
  
  return(list(
    plot = p,
    data = df
  ))
}

###########################################################################



###########################################################################
## Function for donut/stacked pie chart plot
geneStatsPieChart <- function(df, namePrefix, title, fdr_col, lfc_col, fdr_cut = 0.05, lfc_cut = 1){
  
  up_cut <- lfc_cut
  down_cut <- lfc_cut * -1
  
  
  df$significant <- "Non-significant"
  #q value significant but no large fold change
  df[which(df[fdr_col] < p_cutoff), "significant"] <- "Significant"
  
  #large fold change but not low enough q value
  df$upDown <- "non-DEG"
  
  df[which(df[lfc_col] >= up_cut), "upDown"] <- "Up"
  df[which(df[lfc_col] <= down_cut), "upDown"] <- "Down"
  
  
  stats <- df %>% group_by(significant, upDown) %>%
    summarise(n = n())
  
  
  stats <- stats[order(stats$n), ]
  stats$ymax <- cumsum(stats$n)
  stats$ymin <- c(0, head(stats$n, n=-1))
  
  
  p <- ggplot(data = stats) +
    geom_rect(mapping = aes(xmin = 2, xmax = 4, ymin = ymin, ymax = ymax, fill = upDown)) +
    geom_rect(mapping = aes(xmin = 2, xmax = 3, ymin = ymin, ymax = ymax, fill = significant)) +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "non-DEG" = "gray", "Non-significant" = "yellow3", "Significant" = "green")) +
    # coord_cartesian(xlim = c(0,4), ylim = c(0,4000)) + 
    coord_polar("y") +
    theme_bw() +
    xlim(0,4) +
    theme(axis.ticks=element_blank(), axis.text=element_blank(), panel.grid=element_blank())
  
  
  imageName <- paste0(outFilePrefix, "_stats.png" , collapse = "")
  
  png(filename = imageName, width=15000, height=15000, res = 1000)
  print(p)
  dev.off()
  
}

###########################################################################


get_diff_info <- function(degInfoFile, dataPath){
  
  diffInfo <- suppressMessages(readr::read_tsv(degInfoFile))
  
  diffInfo <- dplyr::mutate(
    diffInfo,
    rld = paste(dataPath, "/", comparison, "/", comparison, ".rlogCounts.tab", sep = ""),
    normCount = paste(dataPath, "/", comparison, "/", comparison, ".normCounts.tab", sep = ""),
    deseq2 = paste(dataPath, "/", comparison, "/", comparison, ".DESeq2.tab", sep = ""),
    deseq2_shrink = paste(dataPath, "/", comparison, "/", comparison, ".DESeq2_shrunken.tab", sep = ""),
    deg = paste(dataPath, "/", comparison, "/", comparison, ".DEG_all.txt", sep = ""),
    topGO = paste(dataPath, "/", comparison, "/", comparison, ".topGO.tab", sep = ""),
    keggProfile = paste(dataPath, "/", comparison, "/", comparison, ".keggProfile.tab", sep = ""),
    clusterProfiler_GO = paste(dataPath, "/", comparison, "/", comparison, ".clusterProfiler.GO.tab", sep = ""),
    clusterProfiler_kegg = paste(dataPath, "/", comparison, "/", comparison, ".clusterProfiler.kegg.tab", sep = "")
  )
  
  return(diffInfo)
}








###########################################################################


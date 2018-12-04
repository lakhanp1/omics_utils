




###########################################################################
## function to plot volcano plot
volcanoPlot = function(df, namePrefix, title,
                       fdr_col, lfc_col, fdr_Cut = 0.05, lfc_cut = 1,
                       geneOfInterest, ylimit = 150){
  
  if (is.null(df[["geneName"]])) {
    df$geneName = df$geneId
  }
  
  df$significant = "Non-significant"
  #q value significant but no large fold change
  df[which(df[[fdr_col]] < fdr_Cut & abs(df[[lfc_col]]) < lfc_cut), "significant"] = "Significant"
  
  #large fold change but not low enough q value
  df[which(df[[fdr_col]] > fdr_Cut & abs(df[[lfc_col]]) > lfc_cut), "significant"] = "Fold Change"
  
  #q value significant and large fold change
  df[which(df[[fdr_col]] < fdr_Cut & abs(df[[lfc_col]]) > lfc_cut), "significant"] = "Significant & Fold Change"
  
  #Drwa Volcano plot
  p = ggplot(data = df, mapping = aes_string(x = lfc_col, y = sprintf("-log10(%s)", fdr_col) ) ) +
    geom_point(mapping = aes(color = significant), alpha=0.4, size=1.75) +
    scale_color_manual(values = c("Non-significant" = "black", "Significant" = "red", "Fold Change" = "orange", "Significant & Fold Change" = "green"), name = "Significance") +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, ylimit), expand = TRUE) +
    xlab("log2(fold_change)") + ylab("-log10(q-value)") +
    theme_bw() +
    theme(legend.background = element_rect(colour = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(face = "bold", size = 16),
          # legend.position = c(0.85, 0.85),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          plot.margin = unit(c(1,1,1,1),"cm"),
          axis.title = element_text(size = 28, face = "bold"),
          axis.text = element_text(size = 26)) +
    ggtitle(title)
  
  
  ## highlight genes of interest
  if(!missing(geneOfInterest)){
    tmpDf = dplyr::filter(df, geneName %in% geneOfInterest)
   
    p = p + geom_text_repel(data = tmpDf, mapping = aes(label = geneName),
                            segment.color = '#cccccc',
                            segment.size = 1,
                            size = 5) +
      geom_point(data = tmpDf,
                 color = "black",
                 shape = 1
      )
  }

  
  
  # imageName = paste0(paste0(c(namePrefix, "_volcanoPlot"), collapse = ""), ".png", collapse = "");
  # 
  # png(filename = imageName, width=6000, height=7000, res = 520)
  # print(p)
  # dev.off()
  
  return(p)
}

# plotTitle = paste("Volcano plot for", compare[1], "vs", compare[2], "gene expression comparison", sep = " ")
# p <- volcanoPlot(df = diffData, namePrefix = outFilePrefix, title = plotTitle, fdr_col = "padj", fdr_col = "log2FoldChange")

###########################################################################



###########################################################################
## Function for donut/stacked pie chart plot
geneStatsPieChart <- function(df, namePrefix, title, fdr_col, lfc_col){
  
  p_cutoff <- 0.05
  lfc_cut <- 1
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
  
  
  stats = stats[order(stats$n), ]
  stats$ymax = cumsum(stats$n)
  stats$ymin = c(0, head(stats$n, n=-1))
  
  
  p <- ggplot(data = stats) +
    geom_rect(mapping = aes(xmin = 2, xmax = 4, ymin = ymin, ymax = ymax, fill = upDown)) +
    geom_rect(mapping = aes(xmin = 2, xmax = 3, ymin = ymin, ymax = ymax, fill = significant)) +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "non-DEG" = "gray", "Non-significant" = "yellow3", "Significant" = "green")) +
    # coord_cartesian(xlim = c(0,4), ylim = c(0,4000)) + 
    coord_polar("y") +
    theme_bw() +
    xlim(0,4) +
    theme(axis.ticks=element_blank(), axis.text=element_blank(), panel.grid=element_blank())
  
  
  imageName = paste0(outFilePrefix, "_stats.png" , collapse = "")
  
  png(filename = imageName, width=15000, height=15000, res = 1000)
  print(p)
  dev.off()
  
}

###########################################################################


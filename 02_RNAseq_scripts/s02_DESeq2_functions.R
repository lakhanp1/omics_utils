




###########################################################################
## function to plot volcano plot
#' Generate volcano plot from RNAseq result table
#'
#' @param df data frame
#' @param title plot title
#' @param fdr_col FDR column name
#' @param lfc_col log2(fold change) column name
#' @param fdr_cut FDR cutoff. Default: 0.05
#' @param lfc_cut log2(fold change). Default: 1
#' @param ylimit Y axis limit for volcano plot
#' @param xlimit X axis limit for volcano plot
#' @param markGenes A vector of genes to mark with gene name. Default: NULL
#' @param geneNameCol column from which gene name to use for marking. Default: NULL
#' @param showNames whether to show gene names on volcano plot. Default: TRUE
#' @param highlightGenes A named list of gene IDs for highlighting on volcano plot. Default: NULL
#' @param highlightColor A named vector of colors for each list member. Default: NULL
#'
#' @return A list object with following elements:
#' \itemize{
#' \item plot: A ggplot2 object for volcano plot
#' \item data: Final dataframe used in ggplot2 object
#' }
#' 
#' @export
#'
#' @examples NA
volcano_plot <- function(
  df, title = "volcano plot",
  fdr_col, lfc_col, fdr_cut = 0.05, lfc_cut = 1,
  ylimit = 150, xlimit = c(-4, 4),
  markGenes = NULL, geneNameCol = NULL, showNames = TRUE,
  highlightGenes = NULL, highlightColor = NULL){
  
  if(!is.null(geneNameCol)){
    if (is.null(df[[geneNameCol]])) {
      df[[geneNameCol]] <- df$geneId
    }
  }
  
  up_cut <- lfc_cut
  down_cut <- -1 * lfc_cut
  
  plotDf <- df %>% 
    dplyr::mutate(
      category = dplyr::case_when(
        !! sym(fdr_col) < !! fdr_cut & !! sym(lfc_col) >= !! up_cut ~ "Significant Up",
        !! sym(fdr_col) < !! fdr_cut & !! sym(lfc_col) <= !! down_cut ~ "Significant Down",
        !! sym(fdr_col) < !! fdr_cut ~ "Significant",
        TRUE ~ "Non-significant"
      )
    )
  
  degSummary <- dplyr::group_by(plotDf, category) %>% 
    dplyr::tally() %>% 
    dplyr::mutate(
     str = paste(category, ":", n, sep = " ")
    )
  
  subTitle <- stringr::str_c(degSummary$str, collapse = " || ")
  
  ## squish the value to limits
  plotDf$log10FDR <- -log10(plotDf[[fdr_col]])
  plotDf$log10FDR <- scales::squish(x = plotDf$log10FDR, range = c(0, ylimit))
  plotDf[[lfc_col]] <- scales::squish(x = plotDf[[lfc_col]], range = xlimit)
  
  
  #Drwa Volcano plot
  pt <- ggplot(mapping = aes(x = !! sym(lfc_col), y = log10FDR)) +
    geom_point(
      data = plotDf,
      mapping = aes(color = category), alpha=0.6, size=1.75
    ) +
    scale_color_manual(
      values = c("Significant Up" = "red", "Significant Down" = "green",
                 "Non-significant" = "black", "Significant" = "grey"), 
      name = "Significance")
  
  
  
  ## optionally, color only genes of interest
  if(is.list(highlightGenes)){
    
    if(is.null(highlightColor)){
      stop("No color provided for highlightGenes")
    }
    
    ## prepare the data
    colorGeneDf <- purrr::map_dfr(
      .x = highlightGenes,
      .f = ~ tibble::tibble(geneId = .x),
      .id = "colorGroup"
    ) %>% 
      dplyr::left_join(y = plotDf, by = "geneId")
    
    ## draw the points
    pt <- ggplot(mapping = aes(x = !! sym(lfc_col), y = log10FDR)) +
      geom_point(data = plotDf, color = "grey", alpha = 0.8) +
      geom_point(data = colorGeneDf,
                 mapping = aes(color = colorGroup)) +
      scale_color_manual(values = highlightColor)
    
  }
  
  
  ## mark genes of interest with gene names
  if(!missing(markGenes) && !is.null(markGenes)){
    tmpDf <- dplyr::filter(plotDf, !!sym(geneNameCol) %in% markGenes)
    
    ## draw the points
    pt <- pt +
      geom_point(
        data = tmpDf,
        color = "black",
        size=1.75, shape = 1
      )
    
    ## show the gene lables
    if(isTRUE(showNames)){
      pt <- pt +
        geom_text_repel(
          data = tmpDf,
          mapping = aes(label = !!sym(geneNameCol)),
          segment.color = 'black',
          segment.size = 1,
          ylim = c(-log10(fdr_cut), NA),
          size = 5) 
    }
  }
  
  
  ## theme and plot annotations
  pt <- pt +
    geom_hline(yintercept = -log10(fdr_cut), color = "black", linetype = "dashed") +
    geom_vline(xintercept = -lfc_cut, color = "black", linetype = "dashed") +
    geom_vline(xintercept = lfc_cut, color = "black", linetype = "dashed") +
    scale_x_continuous(name = "log2(fold change)", limits = xlimit, expand = expansion(mult = 0.02)) +
    scale_y_continuous(name = "-log10(q-value)", limits = c(0, ylimit), expand = expansion(mult = 0.02)) +
    guides(color = guide_legend(nrow = 2, byrow = T)) +
    labs(
      title = stringr::str_wrap(title),
      subtitle = subTitle
    ) +
    theme_bw() +
    theme(
      legend.background = element_rect(colour = "black"),
      legend.title = element_text(face = "bold", size = 16),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 13),
      axis.title = element_text(face = "bold", size = 15),
      legend.text = element_text(size = 15),
      legend.key.size = unit(1.2,"cm"),
      plot.margin = unit(c(1,1,1,1),"cm"),
    )
  
  
  return(list(
    plot = pt,
    data = plotDf
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
    fpkm = paste(dataPath, "/", comparison, "/", comparison, ".FPKM.tab", sep = ""),
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


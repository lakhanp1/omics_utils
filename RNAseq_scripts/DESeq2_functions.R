




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
#' @param highlightGenesets A named list of gene IDs for highlighting on volcano plot. Default: NULL
#' @param genesetColor A named vector of colors for each list member. Default: NULL
#' @param pointSize point size in geom_point()
#' @param pointAlpha point alpha in geom_point()
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
  pointSize = 1, pointAlpha = 1,
  markGenes = NULL, geneNameCol = NULL,
  highlightGenesets = NULL, genesetColor = NULL){
  
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
  
  if(!missing(markGenes) && !is.null(markGenes)){
    plotDf <- dplyr::left_join(
      x = plotDf, y = tibble(geneId = markGenes, markGene = TRUE), by = "geneId"
    ) %>% 
      dplyr::mutate(
        geneLabel = dplyr::if_else(
          condition = markGene, true = !!sym(geneNameCol), false = "", missing = ""
        )
      )
  }
  
  # Draw Volcano plot
  pt_base <- ggplot(
    data = plotDf,
    mapping = aes(x = !! sym(lfc_col), y = log10FDR, text = geneId)
  )
  
  
  pt_volc <- pt_base +
    geom_point(
      mapping = aes(color = category), alpha = pointAlpha, size = pointSize, shape = 19
    ) +
    scale_color_manual(
      values = c("Significant Up" = "red", "Significant Down" = "green",
                 "Non-significant" = "black", "Significant" = "grey"), 
      name = "Significance")
  
  ## optionally, color only genes of interest
  if(is.list(highlightGenesets)){
    
    if(is.null(genesetColor)){
      stop("No color provided for highlightGenesets")
    }
    
    ## prepare the data
    colorGeneDf <- purrr::map_dfr(
      .x = highlightGenesets,
      .f = ~ tibble::tibble(geneId = .x),
      .id = "colorGroup"
    ) %>% 
      dplyr::left_join(y = plotDf, by = "geneId")
    
    ## draw the points
    pt_volc <- pt_base +
      geom_point(color = "grey", alpha = pointAlpha, size = pointSize, shape = 19) +
      geom_point(
        data = colorGeneDf,
        mapping = aes(color = colorGroup),
        alpha = pointAlpha, size = pointSize, shape = 19
      ) +
      scale_color_manual(values = genesetColor)
    
  }
  
  
  ## mark genes of interest with gene names
  if(!missing(markGenes) && !is.null(markGenes)){
    
    pt_volc <- pt_volc +
      geom_text_repel(
        mapping = aes(label = geneLabel),
        segment.color = 'black',
        segment.size = 1, min.segment.length = 1,
        xlim  = c(down_cut, up_cut),
        max.overlaps = 5000,
        ylim = c(-log10(fdr_cut), NA),
        size = 5) 
    
    tmpDf <- dplyr::filter(plotDf, geneId %in% markGenes)
    
    ## draw the points
    pt_volc <- pt_volc +
      geom_point(
        data = tmpDf,
        color = "red",
        size = pointSize + 0.5, shape = 1, stroke = 1.5
      )
    
  }
  
  
  ## theme and plot annotations
  pt_volc <- pt_volc +
    geom_hline(yintercept = -log10(fdr_cut), color = "black", linetype = "dashed") +
    geom_vline(xintercept = -lfc_cut, color = "black", linetype = "dashed") +
    geom_vline(xintercept = lfc_cut, color = "black", linetype = "dashed") +
    scale_x_continuous(name = "log2(fold change)", limits = xlimit, expand = expansion(mult = 0.02)) +
    scale_y_continuous(name = "-log10(q-value)", limits = c(0, ylimit), expand = expansion(mult = 0.02)) +
    guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 5))) +
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
      axis.text = element_text(size = 20),
      axis.title = element_text(face = "bold", size = 20),
      panel.grid = element_blank(),
      legend.text = element_text(size = 15),
      plot.margin = unit(c(1,1,1,1),"cm"),
    )
  
  
  return(list(
    plot = pt_volc,
    data = plotDf
  ))
}

###########################################################################



###########################################################################
## Function for donut/stacked pie chart plot
gene_stats_pie_chart <- function(df, namePrefix, title, fdr_col, lfc_col, fdr_cut = 0.05, lfc_cut = 1){
  
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
    kegg = paste(dataPath, "/", comparison, "/", comparison, ".KEGG_fora.tab", sep = ""),
    fgsea = paste(dataPath, "/", comparison, "/", comparison, ".GO_KEGG_fgsea.tab", sep = ""),
    clusterProfiler_GO = paste(dataPath, "/", comparison, "/", comparison, ".clusterProfiler.GO.tab", sep = ""),
    clusterProfiler_kegg = paste(dataPath, "/", comparison, "/", comparison, ".clusterProfiler.kegg.tab", sep = "")
  )
  
  return(diffInfo)
}








###########################################################################


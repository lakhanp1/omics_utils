suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(tm))  # for text mining
suppressPackageStartupMessages(library(SnowballC)) # for text stemming
suppressPackageStartupMessages(library(wordcloud)) # word-cloud generator
suppressPackageStartupMessages(library(configr)) # word-cloud generator
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(KEGGREST))


## This script has functions for using topGO package for GO enrichment
##################################################################################

#' GO enrichment using topGO package
#' 
#' This function performs GO enrichment using Fisher's exact test on gene counts using \code{topGO} package.
#'
#' @param goMapFile gene to GO term assignment file
#' @param genes a vector of geneIds. These geneIds should be present in the first column of goMapFile
#' @param type GO type to use for enrichment. One of BP, MF, CC. Default: BP
#' @param goNodeSize default: 1
#' @param algo One of the algorithm for topGO: "classic", "elim", "weight", "weight01", "lea", "parentchild".
#' Default: "weight01"
#' @param bgNodeLimit Any GO term with more than this number of genes in background is removed from the
#' topGO output. Default: NULL i.e. no such limit is used.
#' @param orgdb org.db for mapping geneIds to gene names
#' @param genenameKeytype gene name column from org.db
#' @param inKeytype org.db keytype for input genes
#' @param topgoKeytype org.db keytype for topGO mapping file geneIds
#'
#' @return A topGO enrichment result table
#' @export
#'
#' @examples topGO_enrichment(goMapFile = goToGeneFile, genes = genes)
#' 
topGO_enrichment <- function(
  goMapFile, genes, type = "BP", goNodeSize = 1, algo = "weight01",
  bgNodeLimit = NULL, orgdb, inKeytype = "GID", topgoKeytype = "GID",
  genenameKeytype){
  
  geneID2GO <- topGO::readMappings(file = goMapFile)
  geneNames <- names(geneID2GO)
  
  
  if(inKeytype != topgoKeytype){
    genes <- suppressMessages(
      AnnotationDbi::mapIds(x = orgdb, keys = genes, column = topgoKeytype, keytype = inKeytype)
    )
  }
  
  genes <- unique(genes)
  
  geneList <- factor(as.integer(geneNames %in% genes))
  names(geneList) <- geneNames
  
  goData <- suppressMessages(
    new(Class = "topGOdata",
        ontology = type,
        allGenes = geneList,
        annot = annFUN.gene2GO,
        gene2GO = geneID2GO,
        nodeSize = goNodeSize)
  )
  
  # nodeCount <- length(attributes(attributes(goData)$graph)$nodes)
  nodeCount <- length(topGO::usedGO(goData))
  
  resFisherWeight <- suppressMessages(
    topGO::runTest(goData, algorithm = algo, statistic = "fisher")
  )
  
  ## get the result table, filter by pValue cutoff 0.05 and calculate Rich factor
  resultTab <- topGO::GenTable(
    goData, 
    pvalue = resFisherWeight,
    orderBy = "resFisherWeight", 
    ranksOf = "resFisherWeight", 
    topNodes = nodeCount, 
    numChar = 1000
  )
  
  
  ## filtering and new column calculations
  resultTab <- resultTab %>%
    dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
    dplyr::filter(pvalue <= 0.05) %>%
    dplyr::rename(
      Observed = Significant
    ) %>% 
    dplyr::mutate(
      richness = as.numeric(sprintf(fmt = "%.3f", (Observed / Annotated)) ),
      log10_pval = as.numeric(sprintf(fmt = "%.3f", -log10(as.numeric(pvalue))))
    ) %>%
    dplyr::select(GO.ID,Term, pvalue, everything())
  
  if(!is.null(bgNodeLimit)){
    resultTab <- dplyr::filter(.data = resultTab, Annotated <= bgNodeLimit)
  }
  
  ## return empty DF if no enrichment
  if(nrow(resultTab) == 0){
    return(resultTab)
  }
  
  ann.genes <- topGO::genesInTerm(object = goData, whichGO = resultTab$GO.ID)
  
  ## prepare mapped gene list
  enrichedGenes <- purrr::map_dfr(
    .x = ann.genes,
    .f = function(x){
      termGenes <- intersect(x, genes)
      return(list(geneIds = paste(termGenes, collapse = ";")))
    },
    .id = "GO.ID"
  )
  
  ## optionally get gene names from orgdb
  if(!missing(orgdb)){
    ## extract gene names
    mappedNames <- suppressMessages(
      AnnotationDbi::select(
        x = orgdb, keys = genes, keytype = topgoKeytype, columns = genenameKeytype)
    ) %>% 
      dplyr::mutate(
        !!sym(genenameKeytype) := if_else(
          condition = is.na(!!sym(genenameKeytype)),
          true = !!sym(topgoKeytype),
          false = !!sym(genenameKeytype)
        )
      ) %>% 
      tibble::deframe()
    
    enrichedGenes <- purrr::map_dfr(
      .x = ann.genes,
      .f = function(x){
        termGenes <- intersect(x, genes)
        return(
          list(
            geneIds = paste(termGenes, collapse = ";"),
            geneNames = paste(mappedNames[termGenes], collapse = ";")
          )
        )
      },
      .id = "GO.ID"
    )
    
  }
  
  
  ## select best GO term from multiple terms which have same genes from input list
  resultTab <- dplyr::left_join(x = resultTab, y = enrichedGenes, by = "GO.ID") %>% 
    dplyr::group_by(geneIds) %>% 
    dplyr::arrange(desc(log10_pval), .by_group = TRUE) %>% 
    dplyr::slice(1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(desc(log10_pval))
  
  resultTab$inputSize <- length(genes)
  
  # dplyr::glimpse(resultTab)
  
  return(resultTab)
}


##################################################################################

#' Bar chart for functional enrichment result
#'
#' @param df a dataframe returned by topGO enrichment function topGO_enrichment
#' @param title title of the plot
#' @param pvalCol pvalue column name used for point color. Default: pvalue
#' @param termCol GO term column name which is used as Y axis Default: Term
#' @param colorCol column name for deciding bar color. Default: category
#' @param countCol column name for gene count. Default: Observed
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples NA
enrichment_bar <- function(df, title, pvalCol = "pvalue", termCol = "Term",
                           colorCol = "category", countCol = "Observed"){
  
  goData <- dplyr::arrange(df, !!sym(colorCol), desc(!!sym(pvalCol)))
  wrap_80 <- wrap_format(80)
  
  goData[[termCol]] <- wrap_80(goData[[termCol]])
  goData[[termCol]] <- sprintf(fmt = "%80s", goData[[termCol]])
  goData[[termCol]] <- factor(goData[[termCol]],
                              levels = unique(goData[[termCol]])
  )
  goData$log10_pval <- -log10(as.numeric(goData[[pvalCol]]))
  
  title <- wrap_80(title)
  
  logPvalCol <- "log10_pval"
  
  ggBar <- ggplot(data = goData, mapping = aes(x = !!sym(termCol), y = !!sym(logPvalCol))) +
    geom_bar(mapping = aes(fill = !!sym(colorCol)), stat='identity') +
    geom_text(mapping = aes(y = !!sym(logPvalCol), label = !!sym(countCol)), hjust = 1.2) +
    labs(title = title, y = "-log10(p-value)") +
    coord_flip() +
    # scale_y_discrete(expand = expand_scale(add = c(0, 1))) +
    facet_grid(rows = vars(!!sym(colorCol)), scales = "free_y") +
    theme_bw() +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.spacing = unit(0, "cm"),
      # panel.border = element_blank(),
      axis.text = element_text(size = 12),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 12, face = "bold", hjust = 1),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  return(ggBar)
}

##################################################################################


#' GO enrichment scatter plot
#' 
#' This function plots the topGO enrichment results in form of a scatter plot. Internally it uses ggplot for plotting.
#'
#' @param df a dataframe returned by topGO enrichment function topGO_enrichment
#' @param title title of the plot
#' @param pvalCol pvalue column name used for point color. Default: pvalue
#' @param termCol GO term column name which is used as Y axis Default: Term
#' @param xVar X axis variable. Default: richness
#' @param sizeVar Point size variable. Default: Observed
#'
#' @return A ggplot object 
#' @export
#'
#' @examples enrichment_scatter(df = goData, title = "topGO enrichment")
#' 
enrichment_scatter <- function(df, title, pvalCol = "pvalue", termCol = "Term",
                               xVar = "richness", sizeVar = "Observed"){
  
  
  goData <- dplyr::arrange(df, desc(!!as.name(pvalCol)))
  wrap_80 <- wrap_format(80)
  goData[[termCol]] <- wrap_80(goData[[termCol]])
  goData[[termCol]] <- sprintf(fmt = "%80s", goData[[termCol]])
  goData[[termCol]] <- factor(goData[[termCol]],
                              levels = unique(goData[[termCol]])
  )
  goData$log10_pval <- -log10(as.numeric(goData[[pvalCol]]))
  
  title <- wrap_80(title)
  logPvalCol <- "log10_pval"
  
  ## color scales
  # scaleLim <- ceiling(min(5, max(goData[[logPvalCol]])))
  scaleLim <- 5
  brk = c(1.30103, 2:scaleLim, scaleLim+1)
  scaleLabels <- c(format(1/(10^c(1.30103, 2:scaleLim)), drop0trailing = T, scientific = F), "smaller")
  
  # structure(l10Vals, names = format(1/(10^seq(1, 3, by = 0.1)), drop0trailing = T, scientific = F))
  
  ## size parameters
  sizeBreaks <- ceiling(seq(from = 1, to = max( max( goData[[sizeVar]] ) + 1, 5), length.out = 5) )
  
  ## ggplot object
  goScatter <- ggplot(data = goData) +
    geom_point(mapping = aes(x = !! as.name(xVar),
                             y = !! as.name(termCol), 
                             size = !! as.name(sizeVar),
                             color = !! as.name(logPvalCol))) +
    scale_color_gradientn(name = "p-value",
                          values = rescale(c(1, 2.5, scaleLim, max(goData[[logPvalCol]], scaleLim+0.1))),
                          colours = c("green", "red", "blue", "darkblue"),
                          breaks = brk,
                          labels = scaleLabels,
                          guide = guide_colorbar(barheight = 10, draw.llim = FALSE, order = 1),
                          oob = squish,
                          limits = c(1, scaleLim + 1)
    ) +
    scale_size_continuous(limits = c(0, max( goData[[sizeVar]] )),
                          range = c(1, 15)) +
    labs(title = title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.8, size = 16, face = "bold"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"))
  
  return(goScatter)
  
}

##################################################################################


#' Perform GO enrichment using topGO and scatter plot
#'
#' @param genes a vector of geneIds. These geneIds should be present in the first column of goMapFile
#' @param goToGeneFile gene to GO term assignment file
#' @param goTitle title for the scatter plot
#' @param plotOut output file name for the scatter plot
#' @param ... Other arguments for topGO_enrichment function
#'
#' @return a topGO results dataframe
#' @export
#'
#' @examples NA
go_and_scatterPlot <- function(genes, goToGeneFile, goTitle, plotOut, ...){
  
  goData <- topGO_enrichment(goMapFile = goToGeneFile, genes = genes, ...)
  topGoScatter <- enrichment_scatter(df = goData, title = goTitle)
  
  # draw Heatmap and add the annotation name decoration
  ht <- max(nrow(goData) * 80, 1500)
  # wd <- (min(max(nchar(as.character(goData$Term))), 80) * 30) * 1.5
  wd <- 3600
  rs <- max(min(wd, ht) / 12, 200)
  
  png(filename = plotOut, width = wd, height = ht, res = rs)
  print(topGoScatter)
  dev.off()
  
  return(goData)
}


##################################################################################


## function to generate the plot and write the matrix. can also be called inside dplyr::do()
#' Title
#'
#' @param genes a vector of geneIds. These geneIds should be present in the first column of goMapFile
#' @param title title for the scatter plot
#' @param outPrefix output file name prefix the scatter plot and topGO result table
#' @param mapFile gene to GO term assignment file
#' @param ... Other arguments for topGO_enrichment function
#'
#' @return a dataframe with image dimentions and other attributes 
#' @export
#'
#' @examples NA
topGO_and_plot_asDf <- function(genes, title, outPrefix, mapFile, ...){
  
  goData <- topGO_enrichment(goMapFile = mapFile, genes = genes, ...)
  
  if(nrow(goData) == 0){
    return(data.frame(height = NA, width = NA, title = NA, res = NA, png = NA, stringsAsFactors = F))
  }
  
  topGoScatter <- enrichment_scatter(df = goData, title = title)
  
  ht <- max(nrow(goData) * 80, 1500)
  wd <- (min(max(nchar(as.character(goData$Term))), 80) * 30) * 1.5
  rs <- max(min(wd, ht) / 12, 200)
  
  pngFile <- paste(outPrefix, "_topGO.png", sep = "")
  tabFile <- paste(outPrefix, "_topGO.tab", sep = "")
  
  png(filename = pngFile, width = wd, height = ht, res = rs)
  print(topGoScatter)
  dev.off()
  
  fwrite(x = goData, file = tabFile, sep = "\t", col.names = T, quote = F)
  
  return(data.frame(height = ht, 
                    width = wd, 
                    res = rs, 
                    count = nrow(goData), 
                    title = title, 
                    png = pngFile,
                    stringsAsFactors = F))
}


##################################################################################


#' Group genes into GO terms at specific level
#'
#' @param genes A vector of gene IDs
#' @param orgdb OrgDb
#' @param goLevel GO graph level at which grouping to be performed. Default: 3
#' @param type GO category. One of "MF", "BP", and "CC". Default: "BP"
#' @param ... Other arguments to groupGO function
#'
#' @return A dataframe in with genes grouped into GO terms at specific level
#' @export
#'
#' @examples NA
clusterProfiler_groupGO <- function(genes, orgdb, goLevel = 3, type = "BP", ...){
  
  grpGo <- suppressMessages(
    clusterProfiler::groupGO(gene = genes,
                             OrgDb = orgdb,
                             ont = type,
                             level = goLevel,
                             ...)
  )
  
  df <- dplyr::mutate_if(.tbl = grpGo@result, .predicate = is.factor, .funs = as.character) %>% 
    dplyr::filter(Count > 0)
  
  return(df)
}



##################################################################################



#' Map genes to GO term of interst
#' 
#' This function maps the gene list onto the GO terms of interest. Internally it
#' uses \code{GO.db} and \code{org.db} object for the organism of interest. GOALL key is used
#' to extract the GO term to gene assignment. This is important as using GOALL
#' ensures that a gene is assigned to a GO term and all its parent terms in GOALL
#' key in org.db
#'
#' @param genes A vector of gene Ids
#' @param goTerms GO term ID vector
#' @param orgdb an org.db object
#' @param inKeytype keytype in org.db for the input genes
#'
#' @return A dataframe with GO term to gene mapping statistics
#' @export
#'
#' @examples NA
GO_map <- function(genes, goTerms, orgdb, inKeytype){
  
  # Ontology(goTerms[1])
  # AnnotationDbi::get(goTerms[1], GO.db::GOBPOFFSPRING)
  
  ## build a GO table 
  goTable <- suppressMessages(
    AnnotationDbi::select(
      x = GO.db,
      keys = goTerms,
      columns = c("GOID", "ONTOLOGY", "TERM"),
      keytype = "GOID"))
  
  
  ## get the total number of genes annotated for each ONTOLOGY category
  ontStats <- suppressMessages(
    AnnotationDbi::select(
      x = orgdb,
      keys = AnnotationDbi::keys(x = orgdb, keytype = inKeytype),
      keytype = inKeytype,
      columns = c("GOALL", "ONTOLOGYALL")
    )
  ) %>% 
    dplyr::group_by(ONTOLOGYALL) %>% 
    dplyr::summarise(n = n_distinct(!!sym(inKeytype))) %>% 
    dplyr::filter(!is.na(ONTOLOGYALL))
  
  
  ## IMP: always use GOALL and no GO. GOALL column is build by considering the
  ## graph architecture of GO.db. If a gene is assigned to a particular GO term,
  ## all its parents GO terms are also annotated for that gene. This may not be
  ## true for GO column
  goData <- suppressMessages(
    AnnotationDbi::select(
      x = orgdb,
      keys = goTerms, columns = c(inKeytype),
      keytype = "GOALL"
    )) %>% 
    dplyr::rename(geneId = !!sym(inKeytype))
  
  ## build summary table
  summaryDf <- dplyr::group_by(goData, GOALL) %>% 
    dplyr::summarise(
      count = length(intersect(x = geneId, y = genes)),
      inputSize = n_distinct(genes),
      background = n_distinct(geneId),
      genes = paste(intersect(x = geneId, y = genes), collapse = ";")
    ) %>% 
    dplyr::ungroup()
  
  goMap <- dplyr::left_join(x = goTable, y = summaryDf,by = c("GOID" = "GOALL")) %>%
    dplyr::left_join(y = ontStats, by = c("ONTOLOGY" = "ONTOLOGYALL")) %>% 
    dplyr::mutate(
      backgroundRatio = paste(background, "/", n, sep = ""),
      enrichment = round(count/background, 3)
    ) %>% 
    dplyr::select(GOID, TERM, ONTOLOGY, count, enrichment, inputSize, backgroundRatio, genes)
  
  return(goMap)
}


##################################################################################


#' Map gene list to KEGG pathways
#'
#' @param genes 
#' @param pathways 
#' @param orgdb 
#'
#' @return
#' @export
#'
#' @examples
KEGG_map <- function(genes, pathways, orgdb){
  
  
}


##################################################################################


#' get GO term IDs at specific level for an ontology
#'
#' @param ont Ontology type. One of \code{c("BP", "CC", "MF")}
#' @param level An integer level. GO tree is traversed from level 1 node.
#'
#' @return A vector of GO IDs at level of interest.
#' @export
#'
#' @examples NA
GO_terms_at_level <- function(ont, level){
  
  ont <- match.arg(arg = ont, choices = c("BP", "CC", "MF"))
  
  switch (ont,
          "BP" = {
            rootNode <- "GO:0008150"
            goChild <- GO.db::GOBPCHILDREN
          },
          "MF" = {
            rootNode <- "GO:0003674"
            goChild <- GO.db::GOMFCHILDREN
          },
          "CC" = {
            rootNode <- "GO:0005575"
            goChild <- GO.db::GOCCCHILDREN
          }
  )
  
  # GOBPCHILDREN, GOMFCHILDREN, GOCCCHILDREN are Bimap objects
  # xx <- as.list(goChild)
  # xx[[rootNode]]
  
  levelNodes <- rootNode
  
  if(level > 1){
    for (i in 2:level) {
      
      levelNodes <- AnnotationDbi::mget(x = c(levelNodes), envir = goChild, ifnotfound = NA) %>% 
        unlist() %>% 
        unique()
      
      levelNodes <- levelNodes[!is.na(levelNodes)]
    }
  }
  
  return(levelNodes)
  
}

##################################################################################

#' Build a geneset of GO terms from orgDb object
#'
#' @param orgdb org.db for mapping gene IDs to GO terms
#' @param outKeytype org.db column name from which geneIds to extract for GO terms
#' @param minNodeSize minimun number of genes annotated to GO term
#' 
#' @inheritParams GO_terms_at_level
#' 
#' @return A named list of genesets for various GO terms
#' @export
#'
#' @examples NA
get_go_geneset <- function(orgdb, outKeytype, level = NULL, ont = "BP", minNodeSize = 5,
                           maxNodeSize = 500){
  ont <- match.arg(arg = ont, choices = c("BP", "CC", "MF"))
  
  if(!is.null(level)){
    goIds <- GO_terms_at_level(ont = ont, level = level)
  } else{
    # AnnotationDbi::keys()
    switch (
      ont,
      "BP" = {
        rootNode <- "GO:0008150"
        goOffspring <- GO.db::GOBPOFFSPRING
      },
      "MF" = {
        rootNode <- "GO:0003674"
        goOffspring <- GO.db::GOMFOFFSPRING
      },
      "CC" = {
        rootNode <- "GO:0005575"
        goOffspring <- GO.db::GOCCOFFSPRING
      }
    )
    
    goIds <- as.list(goOffspring)[[rootNode]]
    
  }
  
  geneToGo <- suppressMessages(
    AnnotationDbi::select(
      x = orgdb, keys = unique(goIds), 
      columns = c(outKeytype, "GOALL"), keytype = "GOALL"
    )) %>% 
    dplyr::filter(!is.na(!!sym(outKeytype)))
  
  geneList <- split(x = geneToGo[[outKeytype]], f = geneToGo$GOALL) %>% 
    purrr::discard(.p = ~ length(.x) < minNodeSize) %>% 
    purrr::discard(.p = ~ length(.x) > maxNodeSize)
  
  return(geneList)
  
}

##################################################################################

#' Gene set from KEGG pathway database
#'
#' @param keggOrg KEGG organism code
#' @param orgdb org.db for mapping gene IDs to GO terms. Default: NULL
#' @param keggKeytype KEGG database geneId keytype from org.db
#' @param outKeytype Output column type from org.db
#'
#' @return A list of genesets belonging to different pathways
#' @export
#'
#' @examples
get_kegg_geneset <- function(keggOrg, orgdb = NULL, keggKeytype = NULL, outKeytype = NULL){
  
  pathways <- tibble::enframe(KEGGREST::keggLink("pathway", keggOrg)) %>% 
    dplyr::rename(keggGeneId = name, pathwayId = value) %>% 
    dplyr::mutate(
      keggGeneId = stringr::str_replace(
        string = keggGeneId, pattern = "\\w+:", replacement = ""
      )
    )
  
  if(!is.null(orgdb)){
    idMap <- suppressMessages(
      AnnotationDbi::select(
        x = orgdb, keys = pathways$keggGeneId, column = outKeytype, keytype = keggKeytype
      )
    )
    
    pathways <- dplyr::left_join(
      x = pathways, y = idMap, by = c("keggGeneId" = keggKeytype)
    ) %>% 
      dplyr::filter(!is.na(!!sym(outKeytype))) %>% 
      dplyr::select(keggGeneId = !!sym(outKeytype), pathwayId)
    
  }
  
  pathwayList <- split(pathways$keggGeneId, pathways$pathwayId)
  
  pathDesc <- tibble::enframe(keggList(database = "pathway", organism = keggOrg)) %>% 
    dplyr::rename(pathwayId = name, description = value)
  
  return(
    list(pathwayList = pathwayList, pathDesc = pathDesc)
  )
  
}

##################################################################################

#' Run fgsea steps to get collapsed list of pathways
#'
#' @param pathways Gene set list
#' @param stats named vector
#' @param ... Other arguments to \code{fgsea()}
#'
#' @return data.table object
#' @export
#'
#' @examples NA
fgsea_steps <- function(pathways, stats, ...){
  # run GSEA analysis
  fgseaRes <- fgsea::fgsea(pathways = pathways, stats = stats, ...)
  
  collapsedPathways <- fgsea::collapsePathways(
    fgseaRes = fgseaRes[order(pval)][pval <= 0.05],
    pathways = pathways,
    stats = stats
  )
  
  uniqueFgsea <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
  
  return(uniqueFgsea)
}


##################################################################################
#' KEGG overrepresentation using \code{fgsea::fora()}
#'
#' @param genes Genes for which KEGG overrepresentation analysis to perform
#' @param keggOrg KEGG organism code
#' @param orgdb org.db for mapping gene IDs
#' @param inKeytype org.db keytype for the input \code{genes}
#' @param keggKeytype org.db keytype for KEGG geneId column
#' @param genenameKeytype org.db keytype for gene name column
#' @param pvalueCutoff P-value cutoff to filter results. Default: 0.05
#' @param ... Other arguments to \code{fgsea::fora()}
#'
#' @return
#' @export
#'
#' @examples
fgsea_kegg_overrepresentation <- function(
  genes, keggOrg, orgdb, inKeytype, keggKeytype,
  genenameKeytype = NULL, pvalueCutoff = 0.05, ...){
  
  keggSet <- get_kegg_geneset(
    keggOrg = keggOrg, orgdb = orgdb,
    keggKeytype = keggKeytype, outKeytype = inKeytype
  )
  
  keggUniverse <- unique(unlist(keggSet$pathwayList))
  # keggUniverse = AnnotationDbi::keys(x = orgdb, keytype = inKeytype),
  
  genes <- intersect(keggUniverse, genes)
  
  if(length(genes) == 0){
    warning("No matching input genes to the universal geneset")
    return(NULL)
  }
  
  enrichRes <- fgsea::fora(
    pathways = keggSet$pathwayList, genes = genes,
    universe = keggUniverse,
    ...
  )
  
  keggDesc <- dplyr::select(enrichRes, pathway) %>% 
    dplyr::left_join(y = keggSet$pathDesc, by = c("pathway" = "pathwayId"))
  
  resultDf <- dplyr::filter(enrichRes, overlap != 0, pval <= pvalueCutoff) %>% 
    dplyr::left_join(y = keggDesc, by = "pathway") %>% 
    dplyr::select(pathway, description, everything())
  
  if(nrow(resultDf) == 0){
    return(NULL)
  }
  
  if(!is.null(genenameKeytype)){
    resultDf <- dplyr::mutate(
      resultDf,
      geneNames = fgsea::mapIdsList(
        x = orgdb, keys = overlapGenes, column = genenameKeytype, keytype = inKeytype
      ),
      richness = as.numeric(sprintf(fmt = "%.3f", (overlap / size)) ),
      inputSize = length(genes)
    ) 
  }
  
  return(resultDf)
}

##################################################################################
#' Perform GSEA on GO term genesets generated from org.db
#'
#' @param genelist A ranked named vector for GSEA analysis
#' @param orgdb org.db for mapping gene IDs
#' @param inKeytype org.db keytype for input \code{genelist}
#' @param level GO tree level
#' @param ont Either BP, CC or MF. Default: BP
#' @param minNodeSize Minimum node size to include GO terms in enrichment analysis
#' @param keggOrg KEGG organism code
#' @param keggKeytype org.db keytype for KEGG geneId column
#' @param genenameKeytype org.db keytype for gene name column
#' @param pvalueCutoff P-value cutoff to filter results. Default: 0.05
#' @param ... Other arguments \code{fgsea::fgsea} function
#' 
#' @return A data.table with fgsea result
#' @export
#'
#' @examples NA
fgsea_orgdb_GO_KEGG <- function(
  genelist, orgdb, inKeytype, level = NULL, ont = "BP", minNodeSize = 5, keggOrg = NULL,
  maxNodeSize = 500, keggKeytype = NULL, genenameKeytype = NULL, pvalueCutoff = 0.05, ...){
  
  # build GO geneset from org.db
  goGenesets <- get_go_geneset(
    orgdb = orgdb, outKeytype = inKeytype, level = level,
    ont = ont, minNodeSize = minNodeSize,
    maxNodeSize = maxNodeSize
  )
  
  gseaGo <- fgsea_steps(pathways = goGenesets, stats = genelist, ...)
  
  ## build a GO table 
  pathDesc <- suppressMessages(
    AnnotationDbi::select(
      x = GO.db,
      keys = gseaGo$pathway,
      columns = c("GOID", "ONTOLOGY", "TERM"),
      keytype = "GOID")
  ) %>% 
    dplyr::filter(!is.na(TERM)) %>% 
    dplyr::select(pathway = GOID, pathway_type = ONTOLOGY, description = TERM)
  
  gseaKegg <- NULL
  
  if(!is.null(keggOrg)){
    
    keggSet <- get_kegg_geneset(
      keggOrg = keggOrg, orgdb = orgdb, keggKeytype = keggKeytype, outKeytype = inKeytype
    )
    
    gseaKegg <- fgsea_steps(pathways = keggSet$pathwayList, stats = genelist, ...)
    
    keggDesc <- dplyr::select(gseaKegg, pathway) %>% 
      dplyr::left_join(y = keggSet$pathDesc, by = c("pathway" = "pathwayId")) %>% 
      dplyr::mutate(pathway_type = "KEGG")
    
    pathDesc <- dplyr::bind_rows(pathDesc, keggDesc)
    
  }
  
  resultDf <- dplyr::bind_rows(gseaGo, gseaKegg) %>% 
    dplyr::left_join(y = pathDesc, by = "pathway") %>% 
    dplyr::filter(pval <= pvalueCutoff) %>% 
    dplyr::arrange(pval)
  
  if(nrow(resultDf) == 0){
    return(NULL)
  }
  
  if(!is.null(genenameKeytype)){
    resultDf <- dplyr::mutate(
      resultDf,
      leadingEdgeNames = fgsea::mapIdsList(
        x = orgdb, keys = leadingEdge, column = genenameKeytype, keytype = inKeytype
      )
    ) 
  }
  
  return(resultDf)
}

##################################################################################

## 
#' find frequent words in the sentences
#'
#' @param sentences A vector of sentences
#' @param topn how many top frequent words to report. Default: all words
#' @param remove remove specific keywords which are common
#'
#' @return
#' @export
#'
#' @examples NA
frequent_key_words = function(sentences, topn = Inf, remove = NULL){
  ## Load the data as a corpus
  docs = Corpus(VectorSource(sentences))
  
  ## filter the data
  docs = tm_map(docs, content_transformer(tolower))
  docs = tm_map(docs, stripWhitespace)
  docs = tm_map(docs, removePunctuation)
  toSpace = content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs = tm_map(docs, toSpace, ":")
  docs = tm_map(docs, toSpace, "-")
  docs = tm_map(docs, removeWords, stopwords("english"))
  docs <- tm_map(docs, removeNumbers)
  
  if (!is.null(remove)) {
    ## c("hypothetical", "protein", "predicted", "putative", "uncharacterized", "nrrl", "cbs", "atcc")
    docs = tm_map(docs, removeWords, remove)
    
  }
  # docs = tm_map(docs, stemDocument) 
  
  dtm = TermDocumentMatrix(docs)
  m = as.matrix(dtm)
  v = sort(rowSums(m),decreasing=TRUE)
  dt = data.frame(word = names(v),freq=v)
  
  if(nrow(dt) < topn){
    topn = nrow(dt)
  }
  
  return(dt)
}


##################################################################################

read_gaf <- function(file){
  gafCols <- c(
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO", "DB_Reference", "EVIDENCE",
    "WithOrFrom", "Aspect", "Name", "Synonym", "DB_Object_Type", "taxon", "Date", "Assigned_by",
    "Annotation_Extension", "Gene_Product_From_ID")
  
  gaf <- suppressMessages(readr::read_tsv(
    file = file, col_names = gafCols, comment = "!"
  ))
  
}


##################################################################################

#' Load TxDb SQLite database as AnnotationDb style object
#'
#' @param org organism code
#' @param yaml YAML config file path which has sqlite file path for different organism codes
#'
#' @return TxDB object
#' @export
#'
#' @examples NA
get_TxDb_sqlite <- function(org, yaml){
  config <- configr::read.config(file = yaml)
  
  file_sqlite <- list.files(
    path = file.path(config[[org]]$TxDb, "inst/extdata"), pattern = "*.sqlite", full.names = T
  )
  
  if(length(file_sqlite) != 1) stop("Unique sqlite file not found")
  
  # txDb <- AnnotationDbi::loadDb(file = file_sqlite, packageName = basename(config[[org]]$TxDb))
  txDb <- suppressMessages(AnnotationDbi::loadDb(file = file_sqlite))
  
  return(txDb)
}


##################################################################################








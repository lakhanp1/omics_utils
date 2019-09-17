library(topGO)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
library(KEGGprofile)
library(clusterProfiler)
library(tm)  # for text mining
library(SnowballC) # for text stemming
library(wordcloud) # word-cloud generator 

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
#' topGO output. Default: 500. If NULL, no such limit is used.
#'
#' @return A topGO enrichment result table
#' @export
#'
#' @examples topGO_enrichment(goMapFile = goToGeneFile, genes = genes)
#' 
topGO_enrichment <- function(goMapFile, genes, type = "BP", goNodeSize = 1,
                             algo = "weight01", bgNodeLimit = 500){
  
  geneID2GO <- topGO::readMappings(file = goMapFile)
  geneNames <- names(geneID2GO)
  
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
  resultTab <- topGO::GenTable(goData, 
                               weightedFisher = resFisherWeight,
                               orderBy = "resFisherWeight", 
                               ranksOf = "weightedFisher", 
                               topNodes = nodeCount, 
                               numChar = 1000)
  
  
  ## filtering and new column calculations
  resultTab <- resultTab %>%
    dplyr::mutate(weightedFisher = as.numeric(weightedFisher)) %>%
    dplyr::filter(weightedFisher <= 0.05) %>%
    dplyr::mutate(richness = as.numeric(sprintf(fmt = "%.3f", (Significant / Annotated)) ),
                  log10_pval = as.numeric(sprintf(fmt = "%.3f", -log10(as.numeric(weightedFisher))))
    ) %>%
    dplyr::select(GO.ID,Term, weightedFisher, everything())
  
  if(!is.null(bgNodeLimit)){
    resultTab <- dplyr::filter(.data = resultTab, Annotated <= bgNodeLimit)
  }
  
  if(nrow(resultTab) == 0){
    return(resultTab)
  }
  
  
  ## add gene name column
  ann.genes <- topGO::genesInTerm(object = goData, whichGO = resultTab$GO.ID)
  
  
  enrichedGenes <- as.data.frame(
    do.call(rbind,
            lapply(X = ann.genes, FUN = function(x){paste(x[x %in% genes], sep = ",", collapse = ",") })
    ),
    stringsAsFactors = F
  )
  
  enrichedGenes <- enrichedGenes %>% 
    tibble::rownames_to_column(var = "GO_ID") %>% 
    dplyr::rename(genes = V1)
  
  ## select best GO term from multiple terms which have same genes from input list
  resultTab <- dplyr::left_join(x = resultTab, y = enrichedGenes, by = c("GO.ID" = "GO_ID")) %>% 
    dplyr::group_by(genes) %>% 
    dplyr::arrange(desc(log10_pval), .by_group = TRUE) %>% 
    dplyr::slice(1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(desc(log10_pval))
  
  resultTab$inputSize <- length(genes)
  return(resultTab)
}


##################################################################################


#' GO enrichment scatter plot
#' 
#' This function plots the topGO enrichment results in form of a scatter plot. Internally it uses ggplot for plotting.
#'
#' @param df a dataframe returned by topGO enrichment function topGO_enrichment
#' @param title title of the plot
#' @param pvalCol pvalue column name used for point color. Default: weightedFisher
#' @param termCol GO term column name which is used as Y axis Default: Term
#' @param xVar X axis variable. Default: richness
#' @param sizeVar Point size variable. Default: Significant
#'
#' @return A ggplot object 
#' @export
#'
#' @examples enrichment_scatter(df = goData, title = "topGO enrichment")
#' 
enrichment_scatter <- function(df, title, pvalCol = "weightedFisher", termCol = "Term",
                               xVar = "richness", sizeVar = "Significant"){
  
  
  goData <- dplyr::arrange(df, desc(!!as.name(pvalCol)))
  wrap_80 <- wrap_format(80)
  goData[[termCol]] <- wrap_80(goData[[termCol]])
  goData[[termCol]] <- sprintf(fmt = "%80s", goData[[termCol]])
  goData[[termCol]] <- factor(goData[[termCol]],
                              levels = unique(goData[[termCol]])
  )
  goData$log10_pval <- -log10(as.numeric(goData[[pvalCol]]))
  
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
#' @param org OrgDb
#' @param goLevel GO graph level at which grouping to be performed. Default: 3
#' @param type GO category. One of "MF", "BP", and "CC". Default: "BP"
#' @param ... Other arguments to groupGO function
#'
#' @return A dataframe in with genes grouped into GO terms at specific level
#' @export
#'
#' @examples NA
clusterProfiler_groupGO <- function(genes, org, goLevel = 3, type = "BP", ...){
  
  grpGo <- suppressMessages(
    clusterProfiler::groupGO(gene = genes,
                             OrgDb = org,
                             ont = type,
                             level = goLevel,
                             ...)
  )
  
  df <- dplyr::mutate_if(.tbl = grpGo@result, .predicate = is.factor, .funs = as.character) %>% 
    dplyr::filter(Count > 0)
  
  return(df)
}


##################################################################################


#' KEGG pathway enrichment using KEGGprofile::find_enriched_pathway
#'
#' @param genes a vector of gene IDs
#' @param orgdb OrgDB database for the organism. 
#' @param keggIdCol Column name from org.db which has ids matching to KEGG database.
#' Usually NCBI gene ids are used as KEGG gene ids.
#' @param keytype Appropriate keytype for the gene list. This keytype is used to extract
#' the geneID to KEGG_gene_id mappings
#' @param keggOrg Three letter organism code. Please refer to KEGG website for details
#' @param pvalCut p-value cutoff for the enriched pathway. Default: 0.05
#' @param qvalCut adjusted p-value cutoff. Default: 1
#' @param minGenes minimum number of annotated genes for enriched pathways. Default: 1
#' @param ... Other arguments to function KEGGprofile::find_enriched_pathway
#'
#' @return A dataframe with enriched pathways
#' @export
#'
#' @examples NA
keggprofile_enrichment <- function(genes, orgdb, keytype, keggIdCol, keggOrg,
                                   pvalCut = 0.05, qvalCut = 1, minGenes = 1, ...){
  
  ## extract KEGG gene IDs 
  keggIds <- suppressMessages(
    AnnotationDbi::select(x = orgdb, keys = genes, keytype = keytype,
                          columns = c(keytype, keggIdCol))
  ) %>% 
    dplyr::filter(!is.na(!!sym(keggIdCol)))
  
  ## KEGG enrichment
  kp <- suppressMessages(
    KEGGprofile::find_enriched_pathway(gene = keggIds[[keggIdCol]], species = keggOrg,
                                       returned_pvalue = pvalCut, returned_adjpvalue = qvalCut,
                                       download_latest = TRUE, returned_genenumber = minGenes,
                                       ...)
  )
  
  
  ## prepare a mapped gene list for the enriched pathways
  assignedGenes <- sapply(X = kp$detail,
                          FUN = function(x){return(paste(x, collapse = ","))},
                          simplify = TRUE, USE.NAMES = TRUE)
  
  assignedDf <- data.frame("pathway_id" = names(assignedGenes),
                           "genes" = assignedGenes, stringsAsFactors = FALSE)
  
  
  ## need to add a modification to return original gene IDs instead of KEGG gene IDs
  
  ## a final result dataframe
  keggDf <- tibble::rownames_to_column(.data = kp$stastic, var = "pathway_id") %>% 
    dplyr::filter(Pathway_Name != "Metabolic pathways") %>% 
    dplyr::left_join(y = assignedDf, by = c("pathway_id" = "pathway_id")) %>% 
    dplyr::rename(Annotated = Gene_Pathway,
                  Significant = Gene_Found,
                  richness = Percentage) %>% 
    dplyr::mutate(inputSize = nrow(keggIds)) %>% 
    dplyr::arrange(pvalue)
  
  return(keggDf)
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
#' @param goTerm A vector of GO term Ids
#' @param org an org.db object
#'
#' @return A dataframe with GO term to gene mapping statistics
#' @export
#'
#' @examples
GO_map <- function(genes, goTerms, org){
  
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
      x = org,
      keys = AnnotationDbi::keys(x = org, keytype = "GID"),
      keytype = "GID",
      columns = c("GOALL", "ONTOLOGYALL")
    )
  ) %>% 
    dplyr::group_by(ONTOLOGYALL) %>% 
    dplyr::summarise(n = n_distinct(GID)) %>% 
    dplyr::filter(!is.na(ONTOLOGYALL))
  
  
  ## IMP: always use GOALL and no GO. GOALL column is build by considering the
  ## graph architecture of GO.db. If a gene is assigned to a particular GO term,
  ## all its parents GO terms are also annotated for that gene. This may not be
  ## true for GO column
  goData <- suppressMessages(
    AnnotationDbi::select(
      x = org,
      keys = goTerms, columns = c("GID"),
      keytype = "GOALL"
    ))
  
  ## build summary table
  summaryDf <- dplyr::group_by(goData, GOALL) %>% 
    dplyr::summarise(
      count = length(intersect(x = GID, y = genes)),
      inputSize = n_distinct(genes),
      background = n_distinct(GID),
      genes = paste(intersect(x = GID, y = genes), collapse = ";")
    ) %>% 
    dplyr::ungroup()
  
  goMap <- dplyr::left_join(x = goTable, y = summaryDf,by = c("GOID" = "GOALL")) %>%
    dplyr::left_join(y = ontStats, by = c("ONTOLOGY" = "ONTOLOGYALL")) %>% 
    dplyr::mutate(backgroundRatio = paste(background, "/", n, sep = "")) %>% 
    dplyr::select(GOID, TERM, ONTOLOGY, count, inputSize, backgroundRatio, genes)
  
  return(goMap)
}


##################################################################################


#' Map gene list to KEGG pathways
#'
#' @param genes 
#' @param pathways 
#' @param org 
#'
#' @return
#' @export
#'
#' @examples
KEGG_map <- function(genes, pathways, org){
  
  
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
#' @examples
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
      
      levelNodes <- mget(x = c(levelNodes), envir = goChild, ifnotfound = NA) %>% 
        unlist() %>% 
        unique()
      
      levelNodes <- levelNodes[!is.na(levelNodes)]
    }
  }
  
  return(levelNodes)
  
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
#' @examples
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
  
  gaf <- data.table::fread(file = file, sep = "\t", header = F,
                           data.table = T, strip.white = T, stringsAsFactors = F,
                           na.strings = "", col.names = gafCols) %>% 
    tibble::as_tibble()
  
  
}











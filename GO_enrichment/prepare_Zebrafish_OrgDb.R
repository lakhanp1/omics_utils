library(tidyverse)
library(data.table)
library(AnnotationForge)
library(GenomicFeatures)
library(biomaRt)

rm(list = ls())


path = "E:/Chris_UM/Database/Zebrafish/GRCz11/annotation_resources"
setwd(path)

##############################################################################


## gene information table
geneInfo <- suppressMessages(
  readr::read_tsv(file = "GRCz11.Ensembl_97.BioMart.gene_info.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                DESCRIPTION = !!as.name("Gene description"),
                GENE_NAME = !!as.name("Gene name"),
                GENE_TYPE = !!as.name("Gene type")) %>% 
  dplyr::mutate(
    DESCRIPTION = if_else(
      condition = is.na(DESCRIPTION), true = GENE_TYPE, false = DESCRIPTION),
    ENSEMBL = GID
  ) %>% 
  tidyr::replace_na(replace = list(GENE_NAME = "")) %>% 
  dplyr::distinct() %>% 
  as.data.frame()


geneInfo$DESCRIPTION <- gsub(pattern = " \\[Source:.*\\]",
                             replacement = "",
                             x = geneInfo$DESCRIPTION, perl = TRUE)


## ensemble to ncbi table
ncbiData <- suppressMessages(
  readr::read_tsv(file = "GRCz11.Ensembl_97.BioMart.NCBI.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                NCBI_ID = !!as.name("NCBI gene ID")) %>% 
  dplyr::filter(!is.na(NCBI_ID)) %>% 
  dplyr::mutate(NCBI_ID = as.character(NCBI_ID)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()

length(unique(ncbiData$GID))


## GO table
goData <- suppressMessages(
  readr::read_tsv(file = "GRCz11.Ensembl_97.BioMart.GO.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                GO = !!as.name("GO term accession"),
                EVIDENCE = !!as.name("GO term evidence code")) %>% 
  dplyr::filter(!is.na(GO)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()

length(unique(goData$GID))

## KEGG table
keggData <- suppressMessages(
  readr::read_tsv(file = "GRCz11.Ensembl_97.BioMart.KEGG.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                KEGG_EC = !!as.name("KEGG Pathway and Enzyme ID")) %>% 
  dplyr::filter(!is.na(KEGG_EC)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()

## ZFIN data
zfinData <- suppressMessages(
  readr::read_tsv(file = "GRCz11.Ensembl_97.BioMart.ZFIN.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                ZFIN_ID = !!as.name("ZFIN ID"),
                ZFIN_SYMBOL = "ZFIN symbol") %>% 
  dplyr::filter(!is.na(ZFIN_ID)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()



makeOrgPackage(
  geneInfo = geneInfo,
  ncbi = ncbiData,
  go = goData,
  kegg = keggData,
  zfin = zfinData,
  version = "10.04.08",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "9606",
  genus = "Danio",
  species = "rerio.GRCz11.Ensembl97",
  goTable = "go",
  verbose = TRUE)













suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(AnnotationForge))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome))

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")

path <- "E:/Chris_UM/Database/A_fumigatus_A1163/annotation_resources/"
setwd(path)


file_chrFeature <- "FungiDB-50_AfumigatusA1163.gene_attributes.tab"
file_goAssociation <- "FungiDB-50_AfumigatusA1163_GO.gaf"
file_orthologsAf293 <- "FungiDB-50_AFumigatus.A1163.Af292.orthologs.txt"

## ANidulans: 162425
## CAlbicans: 5476
## Aspergillus fumigatus: 746128 
## Aspergillus flavus NRRL3357: 332952
## Aspergillus fumigatus A1163: 451804

##############################################################################
## for AFumigatus A1163

geneInfo <- suppressMessages(readr::read_tsv(file = file_chrFeature)) %>% 
  dplyr::rename(
    GID = "Gene ID",
    DESCRIPTION = "Product Description"
  ) %>% 
  dplyr::mutate(
    GENE_NAME = GID
  ) %>% 
  dplyr::select(GID, GENE_NAME, DESCRIPTION)


orthologs <- suppressMessages(readr::read_tsv(file = file_orthologsAf293)) %>% 
  dplyr::rename(
    GID = "Input Ortholog(s)",
    AF293_ID = "Gene ID"
  ) %>% 
  dplyr::select(GID, AF293_ID) %>% 
  tidyr::separate_rows(GID, sep = ",")



##############################################################################
## GO data

goAssociations <- read_gaf(file = file_goAssociation)

## Gene to GO map for topGO
goData =  dplyr::select(goAssociations, DB_Object_ID, GO, EVIDENCE)


## GO dataframe for OrgDb package
goDf = dplyr::left_join(x = geneInfo, y = goData, by = c("GID" = "DB_Object_ID")) %>% 
  dplyr::select(GID, GO, EVIDENCE) %>% 
  dplyr::filter(!is.na(GO)) %>% 
  dplyr::distinct()


## topGO gene-> GO map
topGoMap = dplyr::select(goDf, GID, GO) %>% 
  dplyr::filter(!is.na(GO)) %>% 
  dplyr::group_by(GID) %>%
  dplyr::mutate(GOs = paste(GO, collapse = ",")) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>% 
  dplyr::select(GID, GOs) %>%
  as.data.frame()

readr::write_tsv(x = topGoMap, file = "geneid2go.AFumigatus_A1163.topGO.map", col_names = F)


##############################################################################
## makeOrgPackage

makeOrgPackage(
  geneInfo = geneInfo,
  orthologs = orthologs,
  go = goDf,
  version = "0.0.50",
  maintainer = "Lakhansing Pardeshi <lakhanp@um.edu.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "451804",
  genus = "Aspergillus",
  species = "Fumigatus.A1163",
  goTable = "go",
  verbose = TRUE)


## install package
install.packages("org.AFumigatus.A1163.eg.db", repos = NULL, type = "source")

library(org.AFumigatus.A1163.eg.db)


##############################################################################
file_gff <- "../A_fumigatus_A1163_features.gff"

afuMetadata <- data.frame(
  name = "Resource URL",
  value = "http://www.aspgd.org/"
)

genomeSize <- suppressMessages(
  readr::read_tsv(file = "../genome.size",
                  col_names = c("chr", "length"))) %>% 
  dplyr::mutate(isCircular = FALSE)

seqInfo <- Seqinfo(
  seqnames = genomeSize$chr,
  seqlengths = genomeSize$length,
  isCircular = genomeSize$isCircular,
  genome = "Aspergillus_fumigatus_A1163")


txdbData <- GenomicFeatures::makeTxDbFromGFF(
  file = file_gff,
  dataSource = "A1163 AspGD GFF",
  organism = "Aspergillus fumigatus",
  metadata = afuMetadata,
  taxonomyId = 451804,
  chrominfo = seqInfo
)

makePackageName(txdbData)

makeTxDbPackage(
  txdb = txdbData,
  version = "0.0.50",
  maintainer = "Lakhansing Pardeshi <lakhanp@um.edu.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  destDir = ".",
  pkgname = "TxDb.Afumigatus.A1163.AspGD.GFF"
)

## install package
install.packages("TxDb.Afumigatus.A1163.AspGD.GFF", repos = NULL, type = "source")


columns(txdbData)
exons(txdbData)
transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"))
genes(txdbData)
fiveUtrs <- fiveUTRsByTranscript(txdbData)

##############################################################################
## create BSgenome package
forgeBSgenomeDataPkg(x = "BSgenome.seed", seqs_srcdir = ".")

## install package
install.packages("BSgenome.Afumigatus.A1163.AspGD", repos = NULL, type = "source")



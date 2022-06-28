suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(AnnotationForge))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome))

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")

path <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/annotation_resources/"
setwd(path)


file_chrFeature <- "A_fumigatus_Af293_version_s03-m05-r12_chromosomal_feature.tab"
file_goAssociation <- "FungiDB-45_AfumigatusAf293_GO.gaf"
file_keggId <- "afm_ncbi-geneid.list"


## ANidulans: 162425
## CAlbicans: 5476
## Aspergillus fumigatus: 746128 
## Aspergillus flavus NRRL3357: 332952
## Aspergillus fumigatus A1163: 451804

##############################################################################
## for AFumigatus Af293

geneInfo <- data.table::fread(
  file = file_chrFeature, sep = "\t", header = F, stringsAsFactors = F, na.strings = "", select = c(1:4, 9, 11),
  quote="", col.names = c("GID", "GENE_NAME", "ALIAS", "TYPE", "ASPGD_ID", "DESCRIPTION")) %>% 
  # dplyr::mutate(ALIAS = if_else(condition = is.na(ALIAS), true = GID, false = paste(GID, ALIAS, sep = "|"))) %>% 
  dplyr::mutate(GENE_NAME = if_else(condition = is.na(GENE_NAME), true = GID, false = GENE_NAME)) %>% 
  tidyr::replace_na(list(DESCRIPTION = "NA"))


## gene id alias
geneToAlias <- dplyr::select(geneInfo, GID, ALIAS) %>% 
  dplyr::filter(!is.na(ALIAS)) %>% 
  dplyr::mutate(newAlias = strsplit(ALIAS, split = "|", fixed = T)) %>% 
  dplyr::select(-ALIAS) %>% 
  tidyr::unnest(cols = c(newAlias)) %>% 
  dplyr::rename(ALIAS = newAlias)

## KEGG and NCBI ids
keggData <- readr::read_tsv(file = file_keggId) %>% 
  dplyr::left_join(y = geneToAlias, by = c("KEGG_ID" = "ALIAS")) %>% 
  dplyr::select(GID, KEGG_ID, NCBI_ID) %>% 
  dplyr::filter(!is.na(GID)) %>% 
  as.data.frame()


geneInfo <- dplyr::select(geneInfo, -ALIAS)



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


readr::write_tsv(x = topGoMap, file = "geneid2go.AFumigatus_Af293.topGO.map")


##############################################################################
## makeOrgPackage

makeOrgPackage(
  geneInfo = geneInfo,
  aliasInfo = geneToAlias,
  kegg = keggData,
  go = goDf,
  version = "03.05.12",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "746128",
  genus = "Aspergillus",
  species = "Fumigatus.Af293",
  goTable = "go",
  verbose = TRUE)


## install package
install.packages("org.AFumigatus.Af293.eg.db", repos = NULL, type = "source")

library(org.AFumigatus.Af293.eg.db)


##############################################################################
file_gff <- "../A_fumigatus_Af293_version_s03-m05-r06_features.gff"

afuMetadata <- data.frame(
  name = "Resource URL",
  value = "http://www.aspgd.org/"
)

genomeSize <- suppressMessages(
  readr::read_tsv(file = "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/genome.size",
                  col_names = c("chr", "length"))) %>% 
  dplyr::mutate(isCircular = FALSE)

seqInfo <- Seqinfo(
  seqnames = genomeSize$chr,
  seqlengths = genomeSize$length,
  isCircular = genomeSize$isCircular,
  genome = "Aspergillus_fumigatus_Af293_s03-m05-r06")


txdbData <- GenomicFeatures::makeTxDbFromGFF(
  file = file_gff,
  dataSource = "Af293 AspGD GFF",
  organism = "Aspergillus fumigatus",
  metadata = afuMetadata,
  taxonomyId = 746128,
  chrominfo = seqInfo
)

makePackageName(txdbData)

makeTxDbPackage(
  txdb = txdbData,
  version = "03.05.06",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  destDir = ".",
  pkgname = "TxDb.Afumigatus.Af293.AspGD.GFF"
)



columns(txdbData)
exons(txdbData)
transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"))
genes(txdbData)
fiveUtrs <- fiveUTRsByTranscript(txdbData)

## install package
install.packages("TxDb.Afumigatus.Af293.AspGD.GFF", repos = NULL, type = "source")

##############################################################################
## create BSgenome package
forgeBSgenomeDataPkg(x = "BSgenome.seed", seqs_srcdir = ".")

install.packages("BSgenome.Afumigatus.Af293.AspGD", repos = NULL, type = "source")


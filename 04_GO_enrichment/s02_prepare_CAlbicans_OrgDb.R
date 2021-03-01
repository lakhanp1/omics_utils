suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(AnnotationForge))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome))


rm(list = ls())

path <- "E:/Chris_UM/Database/C_albicans/SC5314_A21/annotation_resources"
setwd(path)


file_caMappings <- "CAlbicans_mapping.txt"
file_goAssociation <- "FungiDB-45_CalbicansSC5314_GO.gaf"

## ANidulans: 162425
## CAlbicans: 5476
## Aspergillus fumigatus: 746128 
## Aspergillus flavus NRRL3357: 332952
## Aspergillus fumigatus A1163: 451804

##############################################################################

caGenes <- data.table::fread(file = file_caMappings, header = T, 
                            sep = "\t", stringsAsFactors = FALSE, na.strings = "") %>% 
  dplyr::rename(GID = ORF19_ID)

view(dfSummary(caGenes))

chrFeature <- dplyr::select(caGenes, GID, CGD_ID, ASSEMBLY22_ID, GENE_NAME, DESCRIPTION) %>% 
  tidyr::replace_na(list(GENE_NAME = ""))

allelaData <- dplyr::select(caGenes, GID, A22_B_ID, A22_A_REGION, A22_B_REGION, A21_REGION, A21_STRAND) %>% 
  tidyr::replace_na(list(A22_B_ID = "", A22_B_REGION = "", A21_REGION = "", A21_STRAND = ""))


keggData <- dplyr::select(caGenes, GID, KEGG_ID) %>% 
  dplyr::filter(!is.na(KEGG_ID))

ncbiData <- dplyr::select(caGenes, GID, NCBI_ID) %>% 
  dplyr::filter(!is.na(NCBI_ID))

##############################################################################
## GO data

goCols <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO", "Reference", "EVIDENCE", "WithOrFrom", "Aspect", "Name", "Synonym", "DB_Object_Type", "taxon", "Date", "Assigned_by")

goAssociations <- read_gaf(file = file_goAssociation)

## Gene to GO map for topGO
goData =  dplyr::select(goAssociations, DB_Object_ID, GO, EVIDENCE)


## GO dataframe for OrgDb package
goDf <- dplyr::left_join(x = chrFeature, y = goData, by = c("ASSEMBLY22_ID" = "DB_Object_ID")) %>% 
  dplyr::select(GID, GO, EVIDENCE) %>% 
  dplyr::filter(!is.na(GO)) %>% 
  dplyr::distinct()


## topGO gene-> GO map
topGoMap <- dplyr::select(goDf, GID, GO) %>% 
  dplyr::filter(!is.na(GO)) %>% 
  dplyr::group_by(GID) %>%
  dplyr::mutate(GOs = paste(GO, collapse = ",")) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>% 
  dplyr::select(GID, GOs) %>%
  as.data.frame()


fwrite(x = topGoMap,
       file = "geneid2go.CAlbicans.topGO.map",
       col.names = F, row.names = F,
       sep = "\t", eol = "\n")


##############################################################################
## makeOrgPackage
## version is same as the CGD version 
makeOrgPackage(
  geneInfo = chrFeature,
  alleleInfo = allelaData,
  go = goDf,
  kegg = keggData,
  ncbi = ncbiData,
  version = "02.09.10",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "5476",
  genus = "Candida",
  species = "albicans.SC5314",
  goTable = "go",
  verbose = TRUE)



## install package
install.packages("org.Calbicans.SC5314.eg.db", repos = NULL, type = "source")

library(org.Calbicans.SC5314.eg.db)

##############################################################################
file_gff <- "../C_albicans_SC5314_A21_current_features.gff"

metadata <- data.frame(
  name = "Resource URL",
  value = "http://www.candidagenome.org"
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
  dataSource = "SC5314 A21 GFF",
  organism = "Candida albicans",
  metadata = metadata,
  taxonomyId = 5476,
  chrominfo = seqInfo
)

makePackageName(txdbData)

makeTxDbPackage(
  txdb = txdbData,
  version = "0.0.50",
  maintainer = "Lakhansing Pardeshi <lakhanp@um.edu.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  destDir = ".",
  pkgname = "TxDb.Calbicans.SC5314.CGD.GFF"
)


columns(txdbData)
exons(txdbData)
transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"))
genes(txdbData)
fiveUtrs <- fiveUTRsByTranscript(txdbData)

install.packages("TxDb.Calbicans.SC5314.CGD.GFF", repos = NULL, type = "source")

##############################################################################
## create BSgenome package
forgeBSgenomeDataPkg(x = "BSgenome.seed", seqs_srcdir = ".")

## need to correct this package. something is wrong with DESCRIPTION file
## install package
install.packages("BSgenome.CAlbicans.SC5314_A21.CGD", repos = NULL, type = "source")







library(tidyverse)
library(data.table)
library(AnnotationForge)
library(GenomicFeatures)
library(BSgenome)


rm(list = ls())

source("E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/s01_topGO_functions.R")

path <- "E:/Chris_UM/Database/A_Nidulans/annotation_resources/"
setwd(path)


file_chrFeature <- "A_nidulans_FGSC_A4_version_s10-m04-r13_chromosomal_feature.tab"
file_goAssociation <- "FungiDB-45_AnidulansFGSCA4_GO.gaf"
file_SM <- "SM_genes.txt"
file_tf <- "TF_genes.txt"

## ANidulans: 162425
## CAlbicans: 5476
## Aspergillus fumigatus: 746128 
## Aspergillus flavus NRRL3357: 332952
## Aspergillus fumigatus A1163: 451804

##############################################################################

geneInfo <- data.table::fread(file = file_chrFeature,
                              sep = "\t", header = F, stringsAsFactors = F, na.strings = "", select = c(1:4, 9, 11),
                              col.names = c("GID", "GENE_NAME", "ALIAS", "TYPE", "ASPGD_ID", "DESCRIPTION")) %>% 
  dplyr::mutate(GENE_NAME = if_else(condition = is.na(GENE_NAME), true = GID, false = GENE_NAME))


## SM cluster information
smGenes <- suppressMessages(readr::read_tsv(file = file_SM)) %>% 
  dplyr::mutate(GID = geneId,
                SM_GENE = geneId) %>% 
  dplyr::select(GID, SM_GENE, SM_CLUSTER, SM_ID) %>% 
  dplyr::distinct()

## transcription factor information
tfGenes <- suppressMessages(readr::read_tsv(file = file_tf)) %>% 
  dplyr::mutate(TF_GENE = geneId, GID = geneId) %>% 
  dplyr::select(GID, TF_GENE)


## make a table of gene and its synonyms which includes all alias as well as gene itself
aliasTable <- dplyr::select(geneInfo, GID, ALIAS) %>%
  dplyr::filter(!is.na(ALIAS)) %>% 
  dplyr::mutate(ALIAS = strsplit(x = ALIAS, split = "\\|")) %>% 
  tidyr::unnest(cols = c(ALIAS)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()


geneInfo <- dplyr::select(geneInfo, -ALIAS)


## kegg to ncbi map
keggMap <- fread(file = "ani_ncbi-gene.list",
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::mutate(geneId = gsub(pattern = "\\.2", replacement = "", x = KEGG_ID))


## merge the kegg map to aspGD geneIds using gene synonyms
## IMP1: there are some smaller genes which were fragment of bigger gene in older version of annotation
## in annotation update, they are merged with the bigger genes. Select the row with bigger gene (rank == 1)
## Eg.: AN2541, AN3435, AN5540, AN6009
## IMP2: there are some genes which are merged from two genes in previous version of annotations.
## such aspGD gene has two different ANx and ANy synonyms. These both can have keggMapping
## remember to remove such AN genes during filtering (n == 1)
keggDf <- dplyr::bind_rows(
  dplyr::left_join(x = aliasTable, y = keggMap, by = c("ALIAS" = "geneId")),
  dplyr::left_join(x = aliasTable, y = keggMap, by = c("GID" = "geneId"))
) %>% 
  dplyr::filter(!is.na(KEGG_ID)) %>%
  dplyr::select(-ALIAS) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(GID) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n == 1) %>% 
  dplyr::select(GID, KEGG_ID, NCBI_ID) %>% 
  as.data.frame()

# setdiff(keggMap$NCBI_ID, keggDf$NCBI_ID)
# setdiff(keggMap$geneId, keggDf$GID)
# 
# setdiff(keggDf$NCBI_ID, keggMap$NCBI_ID)

geneTable <- dplyr::left_join(x = geneInfo, y = keggDf, by = c("GID" = "GID"))


fwrite(x = geneTable, file = "A_nidulans_FGSC_A4.geneTable.tab", sep = "\t", col.names = T, quote = F)


##############################################################################
## GO data

goAssociations <- read_gaf(file = file_goAssociation)

## Gene to GO map for topGO
goData =  dplyr::select(goAssociations, DB_Object_ID, GO, EVIDENCE)

## GO dataframe for OrgDb package
goDf <- dplyr::left_join(x = geneInfo, y = goData, by = c("GID" = "DB_Object_ID")) %>% 
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
       file = "geneid2go.ANidulans.topGO.map",
       col.names = F, row.names = F,
       sep = "\t", eol = "\n")


##############################################################################
## makeOrgPackage

makeOrgPackage(
  geneInfo = geneInfo,
  aliasInfo = aliasTable,
  smInfo = smGenes,
  tfInfo = tfGenes,
  kegg = keggDf,
  go = goDf,
  version = "10.04.13",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "162425",
  genus = "Aspergillus",
  species = "nidulans.FGSCA4",
  goTable = "go",
  verbose = TRUE)


## install package
install.packages("E:/Chris_UM/Database/A_Nidulans/annotation_resources/org.Anidulans.FGSCA4.eg.db",
                 repos = NULL, type = "source")

library(org.Anidulans.FGSCA4.eg.db)


##############################################################################
file_gff <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_version_s10-m04-r03_features.gff"

metadata <- data.frame(
  name = "Resource URL",
  value = "http://www.aspgd.org/"
)

genomeSize <- suppressMessages(
  readr::read_tsv(file = "E:/Chris_UM/Database/A_Nidulans/genome.size",
                  col_names = c("chr", "length"))) %>% 
  dplyr::mutate(isCircular = FALSE)


seqInfo <- Seqinfo(
  seqnames = genomeSize$chr,
  seqlengths = genomeSize$length,
  isCircular = genomeSize$isCircular,
  genome = "Aspergillus_nidulans_s10-m04-r03")

txdbData <- GenomicFeatures::makeTxDbFromGFF(
  file = file_gff,
  dataSource = "A_nidulans AspGD GFF",
  organism = "Aspergillus nidulans",
  metadata = metadata,
  taxonomyId = 162425,
  chrominfo = seqInfo)

makePackageName(txdbData)

makeTxDbPackage(
  txdb = txdbData,
  version = "03.05.06",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  destDir = ".",
  pkgname = "TxDb.Anidulans.FGSCA4.AspGD.GFF"
)



columns(txdbData)
exons(txdbData)
transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"))
GenomicFeatures::genes(txdbData)
fiveUtrs <- fiveUTRsByTranscript(txdbData)

## install package
install.packages("E:/Chris_UM/Database/A_Nidulans/annotation_resources/TxDb.Anidulans.FGSCA4.AspGD.GFF",
                 repos = NULL, type = "source")


##############################################################################
file_gff2 <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_version_s10-m04-r03_features_tRNA_removed.gff"

metadata <- data.frame(
  name = "Resource URL",
  value = "http://www.aspgd.org/"
)

genomeSize <- suppressMessages(
  readr::read_tsv(file = "E:/Chris_UM/Database/A_Nidulans/genome.size",
                  col_names = c("chr", "length"))) %>% 
  dplyr::mutate(isCircular = FALSE)

seqInfo <- Seqinfo(
  seqnames = genomeSize$chr,
  seqlengths = genomeSize$length,
  isCircular = genomeSize$isCircular,
  genome = "Aspergillus_nidulans_s10-m04-r03")

txdbData2 <- GenomicFeatures::makeTxDbFromGFF(
  file = file_gff2,
  dataSource = "A_nidulans AspGD non-tRNA GFF",
  organism = "Aspergillus nidulans",
  metadata = metadata,
  taxonomyId = 162425,
  chrominfo = seqInfo)

makePackageName(txdbData2)

makeTxDbPackage(
  txdb = txdbData2,
  version = "03.05.06",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  destDir = ".",
  pkgname = "TxDb.Anidulans.tRNA.removed"
)


##############################################################################
## create BSgenome package
forgeBSgenomeDataPkg(x = "BSgenome.seed", seqs_srcdir = ".")

## install package
install.packages("E:/Chris_UM/Database/A_Nidulans/annotation_resources/BSgenome.Anidulans.FGSCA4.AspGD",
                 repos = NULL, type = "source")


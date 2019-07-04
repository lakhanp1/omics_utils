library(tidyverse)
library(data.table)
library(AnnotationForge)
library(GenomicFeatures)
library(BSgenome)

rm(list = ls())


path <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_orgDb/"
setwd(path)


file_chrFeature <- "A_fumigatus_Af293_version_s03-m05-r12_chromosomal_feature.tab"
file_goAssociation <- "E:/Chris_UM/Database/GO_associations/gene_association.aspgd.20190117.tab"
file_keggId <- "afm_ncbi-geneid.list"

## using gene coordinates which include UTRs
file_genes <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r12_genes.bed"

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
  tidyr::unnest() %>% 
  dplyr::rename(ALIAS = newAlias)

## KEGG and NCBI ids
keggData <- readr::read_tsv(file = file_keggId) %>% 
  dplyr::left_join(y = geneToAlias, by = c("KEGG_ID" = "ALIAS")) %>% 
  dplyr::select(GID, KEGG_ID, NCBI_ID) %>% 
  dplyr::filter(!is.na(GID)) %>% 
  as.data.frame()


geneBed <- readr::read_tsv(file = file_genes,
                                  col_names = c("CHR", "START", "END", "NAME", "SCORE", "STRAND"))


geneInfo <- dplyr::left_join(x = geneInfo, y = geneBed, by = c("GID" = "NAME")) %>% 
  dplyr::select(GID, GENE_NAME, TYPE, ASPGD_ID, DESCRIPTION, CHR, START, END, STRAND)



##############################################################################
## GO data

goCols = c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO", "Reference", "EVIDENCE", "WithOrFrom", "Aspect", "Name", "Synonym", "DB_Object_Type", "taxon", "Date", "Assigned_by")

goAssociations = data.table::fread(file = file_goAssociation, sep = "\t", header = F, data.table = T, strip.white = T,
                                   stringsAsFactors = F, na.strings = "", col.names = goCols)

## Gene to GO map for topGO
goData = goAssociations %>%
  dplyr::filter(taxon == "taxon:746128") %>%
  dplyr::select(DB_Object_ID, GO, EVIDENCE)


## GO dataframe for OrgDb package
goDf = dplyr::left_join(x = geneInfo, y = goData, by = c("ASPGD_ID" = "DB_Object_ID")) %>% 
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


fwrite(x = topGoMap,
       file = "geneid2go.AFumigatus_Af293.topGO.map",
       col.names = F, row.names = F,
       sep = "\t", eol = "\n")


##############################################################################
## makeOrgPackage

makeOrgPackage(
  geneInfo = geneInfo,
  kegg = keggData,
  go = goDf,
  version = "03.05.12",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "746128",
  genus = "Aspergillus",
  species = "Fumigatus293",
  goTable = "go",
  verbose = TRUE)


## install package
install.packages("E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_orgDb/org.AFumigatus293.eg.db", repos = NULL, type = "source")

library(org.AFumigatus293.eg.db)




##############################################################################
file_gff <- "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r06_features.gff"

afuMetadata <- data.frame(
  name = "Resource URL",
  value = "http://www.aspgd.org/"
)

txdbData <- GenomicFeatures::makeTxDbFromGFF(file = file_gff,
                                             dataSource = "Af293 AspGD GFF",
                                             organism = "Aspergillus fumigatus",
                                             metadata = afuMetadata,
                                             taxonomyId = 746128)

makePackageName(txdbData)

makeTxDbPackage(txdb = txdbData,
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

##############################################################################
## create BSgenome package
forgeBSgenomeDataPkg(x = "BSgenome.seed", seqs_srcdir = ".")



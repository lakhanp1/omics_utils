library(dplyr)
library(data.table)
library(AnnotationForge)
library(clusterProfiler)
library(summarytools)


rm(list = ls())


path <- "E:/Chris_UM/Database/C_albicans/SC5314_A21/CAlbicansA21_OrgDb"
setwd(path)


file_caMappings <- "CAlbicans_mapping.txt"
file_goAssociation <- "E:/Chris_UM/Database/GO_associations/gene_association.cgd"

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

goAssociations <- fread(file = file_goAssociation, sep = "\t", header = F, data.table = T, stringsAsFactors = F,
                       na.strings = "", col.names = goCols)

## Gene to GO map for topGO
goData <- goAssociations %>%
  dplyr::filter(taxon == "taxon:5476") %>%
  dplyr::select(DB_Object_ID, GO, EVIDENCE)


## GO dataframe for OrgDb package
goDf <- dplyr::left_join(x = chrFeature, y = goData, by = c("CGD_ID" = "DB_Object_ID")) %>% 
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
  species = "albicans",
  goTable = "go",
  verbose = TRUE)



## install package
install.packages("E:/Chris_UM/Database/C_albicans/SC5314_A21/CAlbicansA21_OrgDb/org.Calbicans.eg.db",
                 repos = NULL, type = "source")

library(org.Calbicans.eg.db)



library(dplyr)
library(data.table)
library(AnnotationForge)



rm(list = ls())


path = "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/"
setwd(path)


file_chrFeature = "A_nidulans_FGSC_A4_version_s10-m04-r08_chromosomal_feature.tab"
file_goAssociation = "E:/Chris_UM/Database/GO_associations/gene_association.20180604.aspgd.tab"
file_anGenes = "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

## ANidulans: 162425
## CAlbicans: 5476
## Aspergillus fumigatus: 746128 
## Aspergillus flavus NRRL3357: 332952
## Aspergillus fumigatus A1163: 451804

##############################################################################
## for ANidulans

anGenes = data.table::fread(file = file_anGenes, sep = "\t", header = T, stringsAsFactors = F) %>% 
  dplyr::rename(GID = gene, POSITION = Pos, STRAND = Strand, SM_CLUSTER = SM_cluster)


geneInfo = data.table::fread(file = file_chrFeature,
                   sep = "\t", header = F, stringsAsFactors = F, na.strings = "", select = c(1:4, 9, 11),
                   col.names = c("GID", "GENE_NAME", "ALIAS", "TYPE", "ASPGD_ID", "DESCRIPTION")) %>% 
  dplyr::mutate(ALIAS = if_else(condition = is.na(ALIAS), true = GID, false = paste(GID, ALIAS, sep = "|"))) %>% 
  dplyr::mutate(GENE_NAME = if_else(condition = is.na(GENE_NAME), true = GID, false = GENE_NAME))


chrFeature = dplyr::left_join(x = anGenes, y = geneInfo, by = ("GID" = "GID")) %>% 
  dplyr::select(GID, POSITION, STRAND, GENE_NAME, ALIAS, TYPE, ASPGD_ID, DESCRIPTION)

smGenes = dplyr::filter(anGenes, is_SM_gene == TRUE) %>% 
  dplyr::select(GID, SM_CLUSTER)

tfGenes = dplyr::filter(anGenes, is_TF == TRUE) %>% 
  dplyr::mutate(TF = GID) %>% 
  dplyr::select(GID, TF)


## make a table of gene and its synonyms which includes all alias as well as gene itself
geneToAlias = dplyr::select(chrFeature, GID, ALIAS) %>%
  dplyr::group_by(GID) %>% 
  dplyr::slice(1L) %>% 
  dplyr::do(data.frame(
    GID = .$GID,
    synonyms = unlist(strsplit(x = .$ALIAS, split = "\\|")), stringsAsFactors = F)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(rank = if_else(GID == synonyms, true = 1, false = 2)) %>% 
  dplyr::distinct(GID, synonyms, .keep_all = TRUE)


## kegg to ncbi map
keggMap = fread(file = "ani_ncbi-gene.list",
                sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::mutate(GID = gsub(pattern = "\\.2", replacement = "", x = KEGG_ID))


## merge the kegg map to aspGD geneIds using gene synonyms
## IMP1: there are some smaller genes which were fragment of bigger gene in older version of annotation
## in annotation update, they are merged with the bigger genes. Select the row with bigger gene (rank == 1)
## Eg.: AN2541, AN3435, AN5540, AN6009
## IMP2: there are some genes which are merged from two genes in previous version of annotations.
## such aspGD gene has two different ANx and ANy synonyms. These both can have keggMapping
## remember to remove such AN genes during filtering (n == 1)
keggDf = dplyr::left_join(x = geneToAlias, y = keggMap, by = c("synonyms" = "GID")) %>% 
  dplyr::filter(!is.na(KEGG_ID)) %>% 
  dplyr::group_by(GID) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(rank == 1 | n == 1) %>% 
  dplyr::select(GID, KEGG_ID, NCBI_ID) %>% 
  as.data.frame()


geneTable = dplyr::left_join(x = chrFeature, y = keggDf, by = c("GID" = "GID"))


fwrite(x = geneTable, file = "A_nidulans_FGSC_A4.geneTable.tab", sep = "\t", col.names = T, quote = F)


##############################################################################
## GO data

goCols = c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO", "Reference", "EVIDENCE", "WithOrFrom", "Aspect", "Name", "Synonym", "DB_Object_Type", "taxon", "Date", "Assigned_by")

goAssociations = data.table::fread(file = file_goAssociation, sep = "\t", header = F, data.table = T,
                                   stringsAsFactors = F, na.strings = "", col.names = goCols)

## Gene to GO map for topGO
goData = goAssociations %>%
  dplyr::filter(taxon == "taxon:162425") %>%
  dplyr::select(DB_Object_ID, GO, EVIDENCE)


## GO dataframe for OrgDb package
goDf = dplyr::left_join(x = chrFeature, y = goData, by = c("ASPGD_ID" = "DB_Object_ID")) %>% 
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
       file = "geneid2go.ANidulans.topGO.map",
       col.names = F, row.names = F,
       sep = "\t", eol = "\n")


##############################################################################
## makeOrgPackage

makeOrgPackage(
  geneInfo = chrFeature,
  smInfo = smGenes,
  tfInfo = tfGenes,
  kegg = keggDf,
  go = goDf,
  version = "10.04.08",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "162425",
  genus = "Aspergillus",
  species = "nidulans",
  goTable = "go",
  verbose = TRUE)



# ## install package
# install.packages("E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/org.Anidulans.eg.db", repos = NULL, type = "source")
# 
# library(org.Anidulans.eg.db)
# 





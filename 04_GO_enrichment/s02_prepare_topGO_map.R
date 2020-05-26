library(data.table)
library(dplyr)
library(tibble)


## create Gene to GO map file by reading AspGd/CGD GO data file

rm(list = ls())

path = "E:/Chris_UM/Database/GO_associations"
setwd(path)

file_goAssociation = "gene_association.20180604.aspgd.tab"

rawOutFile = "geneid2go.CAlbicans.raw.tab"
mapOutFile = "geneid2go.CAlbicans.map"


goCols = c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "Reference", "Evidence", "WithOrFrom", "Aspect", "Name", "Synonym", "DB_Object_Type", "taxon", "Date", "Assigned_by")

goAssociations = fread(file = file_goAssociation, sep = "\t", header = F, data.table = T, stringsAsFactors = F,
                       na.strings = "", col.names = goCols)

goSlimDescr = fread(file = "GO_SLIM_BP_description.tab", sep = "\t", header = T, data.table = T, stringsAsFactors = F,
                    na.strings = "")


## ANidulans: 162425
## CAlbicans: 5476
## Aspergillus fumigatus: 746128 
## Aspergillus flavus NRRL3357: 332952
## Aspergillus fumigatus A1163: 451804

##############################################################################
## for AFumigatus
## read the AFumigatus geneId to AspGD id information 

afuGeneInfo = fread(file = "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/A_fumigatus_Af293_version_s03-m05-r09_geneInfo.tab", sep = "\t", header = T, select = c("geneId", "AspGDID"))

df = goAssociations %>%
  dplyr::filter(taxon == "taxon:746128") %>%
  dplyr::select(DB_Object_ID, GO_ID) %>%
  dplyr::group_by(DB_Object_ID) %>%
  # dplyr::mutate(GOs = paste(GO_ID, collapse = ",")) %>%
  dplyr::summarise(GOs = paste(GO_ID, collapse = ",")) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(y = afuGeneInfo, by = c("DB_Object_ID" = "AspGDID")) %>% 
  dplyr::select(geneId, GOs)


fwrite(x = df,
       file = "geneid2go.AFumigatus.20180604.map",
       col.names = F, row.names = F,
       sep = "\t", eol = "\n")






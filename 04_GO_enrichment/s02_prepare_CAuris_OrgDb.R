suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(AnnotationForge))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome))



rm(list = ls())


path <- "E:/Chris_UM/Database/Candida_auris_B8441/CAuris_OrgDb"
setwd(path)

file_genes <- "E:/Chris_UM/Database/Candida_auris_B8441/Candida_auris_B8441.info.tab"
file_blastxSummary <- "E:/Chris_UM/Database/Candida_auris_B8441/nr_blastx/CAuris_blastxSummary.tab"
file_annotations <- "E:/Chris_UM/Database/Candida_auris_B8441/functionalAnnotation/Cand_auri_B8441.v1.annot_att_summary.tab"
file_geneToGo <- "E:/Chris_UM/Database/Candida_auris_B8441/functionalAnnotation/CAuris_gene_to_GO.tab"
file_calbOrtho <- "E:/Chris_UM/Database/Candida_auris_B8441/functionalAnnotation/auris_albicans_orthologs.tab"

##############################################################################

chrFeature <- suppressMessages(readr::read_tsv(file = file_genes)) %>%
  dplyr::rename(GID = geneId, POSITION = position, STRAND = strand, DESCRIPTION = description)

blastxData <- suppressMessages(readr::read_tsv(file = file_blastxSummary)) %>% 
  dplyr::rename(GID = qseqid, BLASTX_SUMMARY = frequentWords) %>%
  dplyr::select(GID, BLASTX_SUMMARY) %>%
  tidyr::replace_na(list(BLASTX_SUMMARY = ""))


ann <- suppressMessages(readr::read_tsv(file = file_annotations))

pfamData <- dplyr::select(ann, gene_id, Pfam) %>%
  tidyr::unnest(Pfam = strsplit(x = Pfam, split = ";", fixed = T)) %>%
  dplyr::filter(grepl(pattern = "^PF", x = Pfam, perl = TRUE)) %>%
  tidyr::separate(col = Pfam, into = c("PFAM_ID", "PFAM_TERM"), sep = "\\|", remove = FALSE) %>% 
  dplyr::select(GID = gene_id, PFAM_ID, PFAM_TERM)


orthologData <- suppressMessages(readr::read_tsv(file = file_calbOrtho)) %>%
  dplyr::filter(hasOrtholog == TRUE) %>%
  dplyr::select(GID = CAuris_id, CGD_ORTHOLOG = CGD_ID, ORTHOLOG_SOURCE = ortho_source)


############################################################################## 
goDf <- suppressMessages(readr::read_tsvfread(file = file_geneToGo)) %>%
  dplyr::rename(GID = geneId, GO = GO_id) %>% 
  dplyr::mutate(EVIDENCE = "IEA") %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::distinct()


## topGO gene-> GO map
topGoMap <- dplyr::select(goDf, GID, GO) %>%
  dplyr::group_by(GID) %>%
  dplyr::mutate(GOs = paste(GO, collapse = ",")) %>%
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>%
  dplyr::select(GID, GOs) %>% 
  as.data.frame()


readr::write_tsv(x = topGoMap, file = "geneid2go.CAuris.topGO.map")

############################################################################## 
## makeOrgPackage version is same as the CGD version
makeOrgPackage(geneInfo = chrFeature,
               blastxInfo = blastxData,
               pfamInfo = pfamData,
               albOrtholog = orthologData,
               go = goDf,
               version = "00.00.2",
               maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
               author = "Lakhansing Pardeshi Chris Lab", outputDir = ".", 
               tax_id = "498019",
               genus = "Candida",
               species = "auris",
               goTable = "go",
               verbose = TRUE)


## install package
install.packages("E:/Chris_UM/Database/Candida_auris_B8441/CAuris_OrgDb/org.Cauris.eg.db", repos = NULL, type = "source")

library(org.Cauris.eg.db)

keytypes(org.Cauris.eg.db)

############################################################################## 











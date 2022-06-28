suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(AnnotationForge))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome))


rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")

path <- "D:/work_lakhan/database/Candida_auris_B8441/annotation_resources"
setwd(path)

file_genes <- "Candida_auris_B8441.info.tab"
file_goAssociation <- "FungiDB-54_CaurisB8441_GO.gaf"
file_ebi_goa <- "4144447.C_auris_B8441.goa"
file_annotations <- "../functionalAnnotation/Cand_auri_B8441.v1.annot_att_summary.tab"
file_blastxSummary <- "../nr_blastx/CAuris_blastxSummary.tab"
file_calbOrtho <- "../functionalAnnotation/auris_albicans_orthologs.tab"

##############################################################################

geneInfo <- suppressMessages(readr::read_tsv(file = file_genes)) %>%
  dplyr::rename(GID = geneId, POSITION = position, STRAND = strand, DESCRIPTION = description)

blastxData <- suppressMessages(readr::read_tsv(file = file_blastxSummary)) %>% 
  dplyr::rename(GID = qseqid, BLASTX_SUMMARY = frequentWords) %>%
  dplyr::select(GID, BLASTX_SUMMARY) %>%
  tidyr::replace_na(list(BLASTX_SUMMARY = ""))


ann <- suppressMessages(readr::read_tsv(file = file_annotations))

pfamData <- dplyr::select(ann, gene_id, Pfam) %>%
  dplyr::mutate(Pfam = strsplit(x = Pfam, split = ";", fixed = T)) %>% 
  tidyr::unnest(Pfam) %>%
  dplyr::filter(grepl(pattern = "^PF", x = Pfam, perl = TRUE)) %>%
  tidyr::separate(col = Pfam, into = c("PFAM_ID", "PFAM_TERM"), sep = "\\|", remove = FALSE) %>% 
  dplyr::select(GID = gene_id, PFAM_ID, PFAM_TERM)


orthologData <- suppressMessages(readr::read_tsv(file = file_calbOrtho)) %>%
  dplyr::filter(hasOrtholog == TRUE) %>%
  dplyr::select(GID = CAuris_id, CGD_ORTHOLOG = CGD_ID, ORTHOLOG_SOURCE = ortho_source)


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


readr::write_tsv(x = topGoMap, file = "geneid2go.CAuris_B8441.topGO.map")

uniprotData <- read_gaf(file = file_ebi_goa) %>% 
  dplyr::filter(DB == "UniProtKB") %>% 
  dplyr::select(GID = DB_Object_Symbol, UNIPROT = DB_Object_ID) %>% 
  dplyr::distinct()
  

############################################################################## 
## makeOrgPackage version is same as the CGD version
makeOrgPackage(geneInfo = geneInfo,
               uniprotInfo = uniprotData,
               blastxInfo = blastxData,
               pfamInfo = pfamData,
               albOrtholog = orthologData,
               go = goDf,
               version = "00.00.2",
               maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
               author = "Lakhansing Pardeshi Chris Lab", outputDir = ".", 
               tax_id = "498019",
               genus = "Candida",
               species = "auris.B8441",
               goTable = "go",
               verbose = TRUE)


## install package
install.packages("org.Cauris.B8441.eg.db", repos = NULL, type = "source")

library(org.Cauris.B8441.eg.db)

keytypes(org.Cauris.B8441.eg.db)


##############################################################################
file_gff <- "../Cand_auris_B8441.gff3"

metadata <- data.frame(
  name = "Resource URL",
  value = "https://fungidb.org/fungidb/app/downloads/release-54/CaurisB8441"
)

genomeSize <- suppressMessages(
  readr::read_tsv(file = "../genome.size",
                  col_names = c("chr", "length"))) %>% 
  dplyr::mutate(isCircular = FALSE)


seqInfo <- Seqinfo(
  seqnames = genomeSize$chr,
  seqlengths = genomeSize$length,
  isCircular = genomeSize$isCircular,
  genome = "Candida_auris_B8441")

txdbData <- GenomicFeatures::makeTxDbFromGFF(
  file = file_gff,
  dataSource = "C_auris GFF",
  organism = "Candida auris",
  metadata = metadata,
  taxonomyId = 498019,
  chrominfo = seqInfo)

makePackageName(txdbData)

makeTxDbPackage(
  txdb = txdbData,
  version = "03.05.06",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  destDir = ".",
  pkgname = "TxDb.Cauris.B8441.GFF"
)



columns(txdbData)
exons(txdbData)
transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"))
GenomicFeatures::genes(txdbData)
fiveUtrs <- fiveUTRsByTranscript(txdbData)

## install package
install.packages("TxDb.Cauris.B8441.GFF", repos = NULL, type = "source", verbose = T)

library(TxDb.Cauris.B8441.GFF)

##############################################################################
## create BSgenome package
forgeBSgenomeDataPkg(x = "BSgenome.seed", seqs_srcdir = ".")

## install package
install.packages("BSgenome.Cauris.B8441", repos = NULL, type = "source")









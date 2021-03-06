suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(AnnotationForge))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(BSgenome))

rm(list = ls())


path = "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/"
setwd(path)

##############################################################################
file_gff <- "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/GRCh38p12.gencode.v30.basic.annotation.sorted.gff3"

metadata <- data.frame(
  name = "Resource URL",
  value = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/"
)

genomeSize <- suppressMessages(
  readr::read_tsv(file = "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/genome.size",
                  col_names = c("chr", "length"))) %>% 
  dplyr::mutate(isCircular = FALSE)

seqInfo <- Seqinfo(seqnames = genomeSize$chr, seqlengths = genomeSize$length,
                   isCircular = genomeSize$isCircular, genome = "GRCh38p12.gencode30")

txdbData <- GenomicFeatures::makeTxDbFromGFF(
  file = file_gff,
  dataSource = "GENCODE.v30",
  organism = "Homo sapiens",
  metadata = metadata,
  taxonomyId = 9606,
  chrominfo = seqInfo
)


makePackageName(txdbData)

makeTxDbPackage(
  txdb = txdbData,
  version = "03.05.06",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  destDir = ".",
  pkgname = "TxDb.Hsapiens.GRCh38p12.gencodev30.basic"
)


keytypes(txdbData)
columns(txdbData)
exonsBTx <- exonsBy(x = txdbData, by = "tx")
tx <- transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"))
tx <- transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"),
                  filter = list(tx_name = c("ENST00000489135.5","ENST00000375127.5")))
genes(txdbData)
fiveUtrs <- fiveUTRsByTranscript(txdbData)

AnnotationDbi::select(x = txdbData, keys = "ENSG00000225931",
                      columns = c("TXID", "TXNAME", "TXTYPE", "EXONID", "EXONNAME", "EXONRANK"), keytype = "GENEID")

df <- AnnotationDbi::select(x = txdbData, keys = keys(txdbData, keytype = "TXID"), keytype = "TXID",
                            columns = c("GENEID", "TXTYPE", "TXNAME")) %>% 
  as.data.frame()


##############################################################################


## gene information table
geneInfo <- suppressMessages(
  readr::read_tsv(file = "GRCh38p12.ensembl_release_96.BioMart.gene_info.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                ENSEMBL_VERSION = !!as.name("Gene stable ID version"),
                DESCRIPTION = !!as.name("Gene description"),
                GENE_NAME = !!as.name("Gene name"),
                GENE_TYPE = !!as.name("Gene type")) %>% 
  dplyr::mutate(
    DESCRIPTION = if_else(
      condition = is.na(DESCRIPTION), true = GENE_TYPE, false = DESCRIPTION),
    ENSEMBL = GID
  ) %>% 
  dplyr::distinct() %>% 
  as.data.frame()

geneInfo$DESCRIPTION <- gsub(pattern = " \\[Source:.*\\]",
                             replacement = "",
                             x = geneInfo$DESCRIPTION, perl = TRUE)

## ensemble to ncbi table
ncbiData <- suppressMessages(
  readr::read_tsv(file = "GRCh38p12.ensembl_release_96.BioMart.NCBI.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                NCBI_ID = !!as.name("NCBI gene ID")) %>% 
  dplyr::filter(!is.na(NCBI_ID)) %>% 
  dplyr::mutate(NCBI_ID = as.character(NCBI_ID)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()

## GO table
goData <- suppressMessages(
  readr::read_tsv(file = "GRCh38p12.ensembl_release_96.BioMart.GO.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                GO = !!as.name("GO term accession"),
                EVIDENCE = !!as.name("GO term evidence code")) %>% 
  dplyr::filter(!is.na(GO)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()

## topGO gene-> GO map
topGoMap <- dplyr::select(goData, GID, GO) %>% 
  dplyr::filter(!is.na(GO)) %>% 
  dplyr::group_by(GID) %>%
  dplyr::mutate(GOs = paste(GO, collapse = ",")) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>% 
  dplyr::select(GID, GOs) %>%
  as.data.frame()


fwrite(x = topGoMap,
       file = "geneid2go.HSapiens.GRCh38p12.topGO.map",
       col.names = F, row.names = F,
       sep = "\t", eol = "\n")

## KEGG table
keggData <- suppressMessages(
  readr::read_tsv(file = "GRCh38p12.ensembl_release_96.BioMart.KEGG.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                KEGG_EC = !!as.name("KEGG Pathway and Enzyme ID")) %>% 
  dplyr::filter(!is.na(KEGG_EC)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()

## HGNC gene symbols
hgncData <- suppressMessages(
  readr::read_tsv(file = "GRCh38p12.ensembl_release_96.HGNC_symbol.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
                HGNC_SYMBOL = !!as.name("HGNC symbol")) %>% 
  dplyr::filter(!is.na(HGNC_SYMBOL)) %>% 
  dplyr::distinct() %>% 
  as.data.frame()



makeOrgPackage(
  geneInfo = geneInfo,
  ncbi = ncbiData,
  go = goData,
  kegg = keggData,
  hgnc = hgncData,
  version = "10.04.08",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "9606",
  genus = "Homo",
  species = "Sapiens.gencodev30",
  goTable = "go",
  verbose = TRUE)


##############################################################################
## parse msigDB XML file to get gene set description
file_msigDb_xml <- "E:/Chris_UM/Database/Human/msigdb_v6.2.xml"

msigDb <- xml2::read_xml(x = file_msigDb_xml)

xml_children(msigDb)

geneSets <- xml_find_all(x = msigDb, xpath = ".//GENESET")

xml_attrs(geneSets[[1]]) %>% enframe() %>% as.data.frame()

# STANDARD_NAME, SYSTEMATIC_NAME, DESCRIPTION_BRIEF

geneSetDf <- tibble(
  STANDARD_NAME = xml_attr(x = geneSets, attr = "STANDARD_NAME"),
  DESCRIPTION_BRIEF = xml_attr(x = geneSets, attr = "DESCRIPTION_BRIEF"),
  CATEGORY_CODE = xml_attr(x = geneSets, attr = "CATEGORY_CODE"),
  SUB_CATEGORY_CODE = xml_attr(x = geneSets, attr = "SUB_CATEGORY_CODE")
)

readr::write_tsv(x = geneSetDf, path = "msigDB_geneset_desc.tab")


##############################################################################




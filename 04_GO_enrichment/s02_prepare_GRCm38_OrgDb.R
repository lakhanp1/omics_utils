library(tidyverse)
library(data.table)
library(AnnotationForge)
library(GenomicFeatures)
library(xml2)

rm(list = ls())


path = "E:/Chris_UM/Database/Mouse/GRCm38.99/annotation_resources/"
setwd(path)

##############################################################################
file_gtf <- "E:/Chris_UM/Database/Mouse/GRCm38.99/Mus_musculus.GRCm38.99.chr.gtf"

metadata <- data.frame(
  name = "Resource URL",
  value = "ftp://ftp.ensembl.org/pub/release-99/tsv/mus_musculus/"
)

genomeSize <- suppressMessages(
  readr::read_tsv(file = "E:/Chris_UM/Database/Mouse/GRCm38.99/genome.size",
                  col_names = c("chr", "length"))) %>% 
  dplyr::mutate(isCircular = FALSE)

seqInfo <- Seqinfo(seqnames = genomeSize$chr, seqlengths = genomeSize$length,
                   isCircular = genomeSize$isCircular, genome = "GRCm38p6.99")

txdbData <- GenomicFeatures::makeTxDbFromGFF(
  file = file_gtf,
  format = "gtf",
  dataSource = "GRCm38.p6 Ensembl release-99",
  organism = "Mus musculus",
  metadata = metadata,
  taxonomyId = 10090,
  chrominfo = seqInfo
)


makePackageName(txdbData)

makeTxDbPackage(txdb = txdbData,
                version = "99.38",
                maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
                author = "Lakhansing Pardeshi Chris Lab",
                destDir = ".",
                pkgname = "TxDb.Mmusculus.GRCm38p6.Ensembl99"
)


keytypes(txdbData)
columns(txdbData)
exonsBTx <- exonsBy(x = txdbData, by = "tx")
tx <- transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"))
tx <- transcripts(txdbData, columns = c("TXID", "TXNAME", "TXTYPE"),
                  filter = list(tx_name = c("ENST00000489135.5","ENST00000375127.5")))
genes(txdbData)
fiveUtrs <- fiveUTRsByTranscript(txdbData)

AnnotationDbi::select(x = txdbData, keys = "ENSMUSG00000118636",
                      columns = c("TXID", "TXNAME", "TXTYPE", "EXONID", "EXONNAME", "EXONRANK"),
                      keytype = "GENEID")

df <- AnnotationDbi::select(x = txdbData, keys = keys(txdbData, keytype = "TXID"), keytype = "TXID",
                            columns = c("GENEID", "TXTYPE", "TXNAME")) %>% 
  as.data.frame()


##############################################################################


## gene information table
geneInfo <- suppressMessages(
  readr::read_tsv(file = "Mus_musculus.GRCm38p6.99.biomart.genes.txt")) %>% 
  dplyr::rename(GID = !!as.name("Gene stable ID"),
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
  readr::read_tsv(file = "Mus_musculus.GRCm38.99.entrez.tsv", col_types = "ccccccccc")) %>% 
  dplyr::rename(GID = !!as.name("gene_stable_id"),
                NCBI = !!as.name("xref")) %>% 
  dplyr::filter(!is.na(NCBI), db_name == "EntrezGene") %>% 
  dplyr::select(GID, NCBI) %>% 
  dplyr::distinct() %>% 
  as.data.frame()

## GO table
goData <- suppressMessages(
  readr::read_tsv(file = "Mus_musculus.GRCm38p6.99.biomart.GO.txt")) %>% 
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
       file = "geneid2go.Mmusculus.GRCm38p6.topGO.map",
       col.names = F, row.names = F,
       sep = "\t", eol = "\n")



makeOrgPackage(
  geneInfo = geneInfo,
  ncbi = ncbiData,
  go = goData,
  version = "99.38",
  maintainer = "Lakhansing Pardeshi <lakhanp@umac.mo>",
  author = "Lakhansing Pardeshi Chris Lab",
  outputDir = ".",
  tax_id = "10090",
  genus = "Mus",
  species = "musculus.GRCm38p6.99",
  goTable = "go",
  verbose = TRUE)


##############################################################################




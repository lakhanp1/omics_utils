library(regioneR)
library(org.AFumigatus293.eg.db)
library(TxDb.Afumigatus.Af293.AspGD.GFF)
library(BSgenome.Afumigatus.AspGD.Af293)
library(Biostrings)


txdb <- TxDb.Afumigatus.Af293.AspGD.GFF
genome <- BSgenome.Afumigatus.AspGD.Af293

cdsGrl <- cdsBy(x = txdb, by = "gene")

cdsGr <- unlist(range(cdsGrl))

cdsGr$geneId <- names(cdsGr)

upstreamRegion <- GenomicRanges::promoters(
  x = cdsGr, upstream = 1000, downstream = 0, use.names = TRUE
)

start(upstreamRegion)[which(start(upstreamRegion) < 0)] <- 1
start(upstreamRegion)[which(start(upstreamRegion) == 1)]

upstreamSeq <- BSgenome::getSeq(x = genome, names = upstreamRegion)

upstreamSeq <- upstreamSeq[which(width(upstreamSeq) == 1000)]

gliGenes <- c("Afu6g09630")

gligRegion <- regioneR::toGRanges(A = "Chr6_A_fumigatus_Af293:2,345,058-2,375,633")

glizOvlp <- findOverlaps(query = upstreamRegion, subject = gligRegion)

glizGenes <- upstreamRegion$geneId[glizOvlp@from]

glizGeneProSeq <- upstreamSeq[which(names(upstreamSeq) %in% glizGenes)]

writeXStringSet(x = glizGeneProSeq, filepath = )


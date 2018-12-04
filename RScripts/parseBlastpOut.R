library(data.table)
library(tibble)
library(dplyr)
library(tm)
library(SnowballC)


rm(list = ls())


path <- "E:/Chris_UM/Database/Candida_auris_B8441/nr_blastx"
setwd(path)


###########################################################################
blastpOutFile = "E:/Chris_UM/Database/Candida_auris_B8441/nr_blastx/Candida_auris_B8441.genes.blastx.tab"
blastpCols = c("qseqid", "qlen", "sseqid", "sgi", "sacc", "pident", "length", "qcovs", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid", "stitle")

bioprojMapFile = "E:/Chris_UM/Database/C_albicans/PRJNA14005/PRJNA14005_prot_to_CGD.tab"
caMapFile = "E:/Chris_UM/Database/C_albicans/SC5314_A21/C_albicans_mappings.txt"
calbicansHSP90DelFile = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/HSP90_nComms_diff/CalbicansHSP90delFoldChange.tab"
degFile = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/noDOX_withDOX_diff_final/noDOX_withDOX_DEG_all.txt"


blastp = fread(file = blastpOutFile, sep = "\t", stringsAsFactors = F, header = F, col.names = blastpCols, na.strings = "N/A")

bioprojMap = fread(file = bioprojMapFile, sep = "\t", stringsAsFactors = F, header = T, drop = c("accession", "locusTag"))
calbicansMap = fread(file = caMapFile, sep = "\t", stringsAsFactors = F, header = T, select = c("CGD_Primary_ID", "CGD_old_id", "ORF19_ID", "GENE_NAME", "Description"))
albicansHSP90Fc = fread(file = calbicansHSP90DelFile, sep = "\t", stringsAsFactors = F, header = T)
degData = fread(file = degFile, sep = "\t", stringsAsFactors = F, header = T)

###########################################################################




###########################################################################
## summarize blastx results

## find frequent words in the description of blastx hits
frequent_key_words = function(sentences, topn){
  ## Load the data as a corpus
  docs <- Corpus(VectorSource(sentences))
  
  ## filter the data
  docs = tm_map(docs, content_transformer(tolower))
  docs = tm_map(docs, stripWhitespace)
  docs = tm_map(docs, removePunctuation)
  toSpace = content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs = tm_map(docs, toSpace, ":")
  docs = tm_map(docs, removeWords, stopwords("english"))
  docs = tm_map(docs, removeWords, c("hypothetical", "protein", "predicted", "putative", "uncharacterized", "nrrl", "cbs", "atcc", "aspergillus", "candida", "albicans", "auris", "saccharomyces", "cerevisiae"))
  docs <- tm_map(docs, removeNumbers)
  # docs = tm_map(docs, stemDocument) 
  
  dtm = TermDocumentMatrix(docs)
  m = as.matrix(dtm)
  v = sort(rowSums(m),decreasing=TRUE)
  d = data.frame(word = names(v),freq=v)
  
  if(nrow(d) < topn){
    topn = nrow(d)
  }
  
  freqWords = paste(d$word[1:topn], collapse = "; ")
  return(freqWords)
}


## summarize blast result for each query
blastxSummary = blastp %>%
  dplyr::group_by(qseqid) %>%
  dplyr::summarise(
    nhits = n(),
    minBitscore = min(bitscore),
    maxBitscore = max(bitscore),
    frequentWords = frequent_key_words(sentences = stitle, topn = 10)
  )


fwrite(x = blastxSummary, file = "CAuris_blastxSummary.tab", sep = "\t", quote = F, col.names = T, na = "NA")

###########################################################################




###########################################################################
## C. auris: 498019
## C. albicans: 237561

## filter the blastp results for species of interest
caBestHits = blastp %>%
  dplyr::filter(staxid == "237561") %>%
  dplyr::left_join(y = bioprojMap, by = c("sgi" = "GI")) %>%
  dplyr::left_join(y = calbicansMap, by = c("cgdId" = "CGD_Primary_ID")) %>%
  dplyr::filter(!is.na(cgdId)) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>%
  dplyr::select(qseqid, bitscore, stitle, cgdId, CGD_old_id, ORF19_ID, GENE_NAME, Description) %>%
  dplyr::rename(ca_bitscore = bitscore, ca_stitle = stitle)


caHitsWithSum = dplyr::left_join(x = blastxSummary, y = caBestHits, by = c("qseqid" = "qseqid")) %>%
  dplyr::select(qseqid, nhits, minBitscore, maxBitscore, frequentWords, everything())
  

fwrite(x = caHitsWithSum, file = "CAlbicansBestHits.tab", sep = "\t", quote = F, col.names = T, na = "NA")


degCaOrth = dplyr::filter(degData, diff != "noDEG") %>%
  dplyr::left_join(y = caHitsWithSum, by = c("geneID" = "qseqid")) %>%
  dplyr::left_join(y = albicansHSP90Fc, by = c("ORF19_ID" = "geneId")) %>%
  dplyr::arrange(desc(diff), desc(log2FoldChange))

fwrite(x = degCaOrth, file = "DEG_CAlbicansBestHits.tab", sep = "\t", quote = F, col.names = T, na = "NA")


###########################################################################



###########################################################################
## Saccharomyces cerevisiae: 559292

scBioprojMapFile = "E:/Chris_UM/Database/Yeast_S288C/PRJNA128/PRJNA128_prot_to_SGD.tab"
scFeatureFile = "E:/Chris_UM/Database/Yeast_S288C/annotation/SGD_features.tab"

scBioprojMap = fread(file = scBioprojMapFile, sep = "\t", stringsAsFactors = F, header = T, drop = c("accession", "locusTag"))

sgdData = fread(file = scFeatureFile, sep = "\t", stringsAsFactors = F, header = T, select = c("SGDID", "geneName", "description"))

## filter the blastp results for species of interest
scBestHits = blastp %>%
  dplyr::filter(staxid == "559292") %>%
  dplyr::left_join(y = scBioprojMap, by = c("sgi" = "GI")) %>%
  dplyr::left_join(y = sgdData, by = c("dbId" = "SGDID")) %>%
  dplyr::filter(!is.na(dbId)) %>%
  dplyr::group_by(qacc) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>%
  dplyr::select(qacc, bitscore, stitle, dbId, geneName, description) %>%
  dplyr::rename(sc_bitscore = bitscore, sc_stitle = stitle, sgdId = dbId, scGene = geneName)


degScOrth = dplyr::filter(degData, diff != "noDEG") %>%
  left_join(y = blastxSummary, by = c("geneID" = "qacc")) %>%
  dplyr::left_join(y = scBestHits, by = c("geneID" = "qacc")) %>%
  dplyr::arrange(desc(diff), desc(log2FoldChange))


fwrite(x = degScOrth, file = "DEG_SCerevisiaeBestHits.tab", sep = "\t", quote = F, col.names = T, na = "NA")

###########################################################################




###########################################################################
## after finalizing orthologs, merge the data together
# caHsp90Data = dplyr::left_join(x = calbicansMap, y = albicansHSP90Fc, by = c("ORF19_ID" = "geneId"))
# 
# fwrite(x = caHsp90Data, file = "caHsp90Data.tab", sep = "\t", quote = F, col.names = T, na = "NA")


## merge 

## filter the DEGs with lower fold_change cutoff of 0.585
cutoff = 0.585

upDegs = dplyr::filter(degData, log2FoldChange > cutoff & padj < 0.05)

upDegs$position





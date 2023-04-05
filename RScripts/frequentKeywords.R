
library(data.table)
library(tibble)
library(dplyr)
library(tm)
library(SnowballC)


rm(list = ls())



## find frequent words in the description of blastx hits
frequent_key_words = function(sentences, topn){
  ## Load the data as a corpus
  docs = Corpus(VectorSource(sentences))
  
  ## filter the data
  docs = tm_map(docs, content_transformer(tolower))
  docs = tm_map(docs, stripWhitespace)
  docs = tm_map(docs, removePunctuation)
  toSpace = content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs = tm_map(docs, toSpace, ":")
  docs = tm_map(docs, toSpace, "-")
  docs = tm_map(docs, removeWords, stopwords("english"))
  docs = tm_map(docs, removeWords, c("hypothetical", "protein", "predicted", "putative", "uncharacterized", "nrrl", "cbs", "atcc"))
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



genes = read.table(file = "clipboard", sep = "\t", header = T, stringsAsFactors = F, strip.white = T)
  
frequentWords = frequent_key_words(sentences = genes$genes, topn = 50)


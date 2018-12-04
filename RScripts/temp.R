
library(rtracklayer)
library(Gviz)


rm(list = ls())


path <-  "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data"
setwd(path)


bw <- "bwFiles/An_cclA_20h_HA_1_normalized.bw"


bwgr <- import(con = bw, format = "bigWig", as = "GRanges")

# dat <- sin(seq(pi, 10*pi, len = 500))
# 
# dTrack.big <- DataTrack(
#   start = seq(1, 1e+05, len = 500),
#   width = 15,
#   chromosome = "chrX",
#   genome = "hg19",
#   name = "sinus",
#   data = sin(seq(pi, 5*pi, len = 500))*runif(500, 0.5, 1.5))
# 
# plotTracks(dTrack.big, type = "hist")
# 
# plotTracks(dTrack.big, type = "hist", window = "fixed", windowSize = 1000)


options(ucscChromosomeNames=FALSE)

dt <- DataTrack(range = bwgr, genome = "ANidulans", name = "An_cclA_20h_HA_1")

chromosome(dt)
seqinfo(dt)
seqlevels(dt)

plotTracks(trackList = dt,
           type = "heatmap",
           window = "fixed", windowSize = 1000,
           gradient = c("black", "yellow","green", "red"),
           # aggregation = function(x){return(quantile(x, 0.99, names = FALSE, na.rm = TRUE))}
           )











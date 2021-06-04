# 
# rm(list = ls())
# 
# ############################################################################
# ###############         uninstall packages          ########################
# ############################################################################
# 
# create a list of all installed packages
# ip <- as.data.frame(installed.packages())
# head(ip)
# 
# # if you use MRO, make sure that no packages in this library will be removed
# ip <- subset(ip, !grepl("MRO", ip$LibPath))
# 
# # we don't want to remove base or recommended packages either\
# ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
# 
# # determine the library where the packages are installed
# path.lib <- unique(ip$LibPath)
# 
# # create a vector with all the names of the packages you want to remove
# pkgs.to.remove <- ip[,1]
# head(pkgs.to.remove)
# 
# # remove the packages
# sapply(pkgs.to.remove, remove.packages, lib = path.lib)


############################################################################
###############         install packages            ########################
############################################################################

rm(list = ls())

chooseCRANmirror(graphics = FALSE, ind = 18)

install.packages(
  pkgs = c(
    "installr", "Rcpp", "tidyverse", "ggpubr", "data.table", "configr", "renv",
    "lazyeval", "dendsort", "dendextend", "dynamicTreeCut", "RColorBrewer",
    "hashmap", "reshape", "FactoMineR", "factoextra", "VennDiagram", "imputeTS",
    "summarytools", "UpSetR", "esquisse", "corrgram", "here", "matrixStats",
    "NbClust", "DT", "msigdbr", "openxlsx", "ggbeeswarm", "PoiClaClu", "tm",
    "wordcloud", "SnowballC", "gginnards"
  ),
  dependencies = T)

#Time series data imputation
install.packages("imputeTS", dependencies = T)

install.packages("BiocManager")
BiocManager::install()

## core bioc packages
BiocManager::install(
  pkgs = c(
    "BiocGenerics", "S4Vectors", "IRanges", "GenomicRanges", "Biostrings",
    "AnnotationForge", "GenomeInfoDbData", "AnnotationHub", "AnnotationDbi",
    "GenomicFeatures", "rtracklayer", "BSgenome", "GenomicAlignments",
    "BiocParallel", "SummarizedExperiment", "plyranges", "Rsamtools"
  )
)

BiocManager::install(
  pkgs = c(
    "DESeq2", "tximport", "topGO", "DiffBind", "regioneR",
    "ballgown", "pathview", "DO.db", "KEGGprofile", "preprocessCore", "GO.db",
    "clusterProfiler", "apeglm", "Gviz", "ggbio", "KEGGREST", "KEGG.db"
  )
)


install.packages(c("devtools"), dependencies = T)

library(remotes)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

remotes::install_github("jokergoo/circlize")
remotes::install_github("jokergoo/ComplexHeatmap")
remotes::install_github("jokergoo/EnrichedHeatmap")

# BiocManager::install("ComplexHeatmap")
# BiocManager::install("EnrichedHeatmap")


## install my own packages
library(devtools)

devtools::document(pkg = "E:/Chris_UM/GitHub/markPeaks")
devtools::install("E:/Chris_UM/GitHub/markPeaks", upgrade = "never")
# remotes::install_github(repo = "lakhanp1/markPeaks", ref = "dev", upgrade = "never")

# devtools::load_all(path = "E:/Chris_UM/GitHub/chipmine", reset = TRUE)
devtools::document(pkg = "E:/Chris_UM/GitHub/chipmine")
devtools::install("E:/Chris_UM/GitHub/chipmine", upgrade = "never")
# devtools::install("E:/Chris_UM/GitHub/chipmine")
# remotes::install_github(repo = "lakhanp1/chipmine", ref = "dev_v2", upgrade = "never")
# remotes::install_github(repo = "lakhanp1/chipmine", ref = "1.6.0", upgrade = "never")

##########################################################################
###############         Org.Db packages            #######################
##########################################################################
## A. nidulans
remotes::install_github(
  "lakhanp1/fungal_resources/A_nidulans/org.Anidulans.FGSCA4.eg.db",
  upgrade = "never")
remotes::install_github(
  "lakhanp1/fungal_resources/A_nidulans/TxDb.Anidulans.FGSCA4.AspGD.GFF",
  upgrade = "never")
remotes::install_github(
  "lakhanp1/fungal_resources/A_nidulans/BSgenome.Anidulans.FGSCA4.AspGD",
  upgrade = "never")
# remotes::install(
#   "E:/Chris_UM/Database/A_Nidulans/annotation_resources/TxDb.Anidulans.tRNA.removed",
#   upgrade = "never")


## C. albicans
remotes::install_github(
  "lakhanp1/fungal_resources/C_albicans/org.Calbicans.SC5314.eg.db",
  upgrade = "never")
remotes::install_github(
  "lakhanp1/fungal_resources/C_albicans/TxDb.Calbicans.SC5314.CGD.GFF",
  upgrade = "never")
remotes::install_github(
  "lakhanp1/fungal_resources/C_albicans/BSgenome.CAlbicans.SC5314_A21.CGD",
  upgrade = "never")


## C. auris
remotes::install_github("lakhanp1/fungal_resources/C_auris/org.Cauris.eg.db",
                         upgrade = "never")


## A. fumigatus Af293
remotes::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_Af293/org.AFumigatus.Af293.eg.db",
  upgrade = "never")
remotes::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_Af293/TxDb.Afumigatus.Af293.AspGD.GFF",
  upgrade = "never")
remotes::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_Af293/BSgenome.Afumigatus.Af293.AspGD",
  upgrade = "never")


## A. fumigatus A1163
remotes::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_A1163/org.AFumigatus.A1163.eg.db",
  upgrade = "never")
remotes::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_A1163/TxDb.Afumigatus.A1163.AspGD.GFF",
  upgrade = "never")
remotes::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_A1163/BSgenome.Afumigatus.A1163.AspGD",
  upgrade = "never")


## Human
devtools::install(
  "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/TxDb.Hsapiens.GRCh38p12.gencodev30.basic",
  upgrade = "never")
devtools::install(
  pkg = "E:/Chris_UM/Database/Human/GRCh38p12.gencode30/annotation_resources/org.HSapiens.gencodev30.eg.db",
  upgrade = "never")


## Zebrafish GRCz11
devtools::install(
  pkg = "E:/Chris_UM/Database/Zebrafish/GRCz11/annotation_resources/org.DRerio.GRCz11.Ensembl97.eg.db",
  upgrade = "never")


## Mouse
devtools::install(
  pkg = "E:/Chris_UM/Database/Mouse/GRCm38.99/annotation_resources/TxDb.GRCm38p6.Ensembl100",
  upgrade = "never")
devtools::install(
  pkg = "E:/Chris_UM/Database/Mouse/GRCm38.99/annotation_resources/org.GRCm38p6.Ensembl100.eg.db",
  upgrade = "never")








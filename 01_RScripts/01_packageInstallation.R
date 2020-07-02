# 
# rm(list = ls())
# 
# ############################################################################
# ###############         uninstall packages          ########################
# ############################################################################
# 
# # create a list of all installed packages
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

chooseCRANmirror(graphics = FALSE, ind = 14)

install.packages(
  pkgs = c("installr", "Rcpp", "tidyverse", "ggpubr", "data.table", "readxl",
           "lazyeval", "dendsort", "dendextend", "dynamicTreeCut", "RColorBrewer",
           "hashmap", "reshape", "FactoMineR", "factoextra", "VennDiagram",
           "imputeTS", "summarytools", "UpSetR", "esquisse", "corrgram", "here",
           "matrixStats", "NbClust", "DT", "msigdbr", "openxlsx", "ggbeeswarm"),
  dependencies = T)


if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install()

## core bioc packages
BiocManager::install(
  pkgs = c("BiocGenerics", "S4Vectors", "IRanges", "GenomicRanges", "Biostrings", "GenomeInfoDbData",
           "AnnotationHub", "AnnotationDbi", "GenomicFeatures", "rtracklayer", "BSgenome",
           "GenomicAlignments", "BiocParallel", "SummarizedExperiment", "plyranges")
)

BiocManager::install(
  pkgs = c("Rsamtools", "DESeq2", "tximport", "topGO", "DiffBind", "regioneR",
           "ballgown", "pathview", "DO.db", "KEGGprofile", "preprocessCore", 
           "clusterProfiler", "apeglm", "Gviz", "KEGGREST", "KEGG.db"
  )
)


#Time series data imputation
install.packages("imputeTS", dependencies = T)

install.packages(c("devtools"), dependencies = T)
library(devtools)
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jokergoo/EnrichedHeatmap")

# BiocManager::install("ComplexHeatmap")
# BiocManager::install("EnrichedHeatmap")


## install my own packages
library(devtools)
# devtools::load_all(path = "E:/Chris_UM/GitHub/chipmine", reset = TRUE)
devtools::document(pkg = "E:/Chris_UM/GitHub/chipmine")
devtools::install("E:/Chris_UM/GitHub/chipmine", upgrade = "never")
# devtools::install("E:/Chris_UM/GitHub/chipmine")
# devtools::install_github(repo = "lakhanp1/chipmine", ref = "dev_v.2")

##########################################################################
###############         Org.Db packages            #######################
##########################################################################
## A. nidulans
devtools::install_github(
  "lakhanp1/fungal_resources/A_nidulans/org.Anidulans.FGSCA4.eg.db",
  upgrade = "never")
devtools::install_github(
  "lakhanp1/fungal_resources/A_nidulans/TxDb.Anidulans.FGSCA4.AspGD.GFF",
  upgrade = "never")
devtools::install_github(
  "lakhanp1/fungal_resources/A_nidulans/BSgenome.Anidulans.FGSCA4.AspGD",
  upgrade = "never")
devtools::install(
  "E:/Chris_UM/Database/A_Nidulans/annotation_resources/TxDb.Anidulans.tRNA.removed",
  upgrade = "never")

## C. albicans
devtools::install_github(
  "lakhanp1/fungal_resources/C_albicans/org.Calbicans.SC5314.eg.db",
  upgrade = "never")

## C. auris
devtools::install_github("lakhanp1/fungal_resources/C_auris/org.Cauris.eg.db",
                         upgrade = "never")

## A. fumigatus
devtools::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_Af293/org.AFumigatus.Af293.eg.db",
  upgrade = "never")
devtools::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_Af293/TxDb.Afumigatus.Af293.AspGD.GFF",
  upgrade = "never")
devtools::install_github(
  "lakhanp1/fungal_resources/A_fumigatus_Af293/BSgenome.Afumigatus.Af293.AspGD",
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
  pkg = "E:/Chris_UM/Database/Mouse/GRCm38.99/annotation_resources/TxDb.Mmusculus.GRCm38p6.Ensembl100",
  upgrade = "never")

devtools::install(
  pkg = "E:/Chris_UM/Database/Mouse/GRCm38.99/annotation_resources/org.Mmusculus.GRCm38p6.Ensembl100.eg.db",
  upgrade = "never")


############################################################################
###############         BSgenome packages            #######################
############################################################################

# ## install A_fumigatus_293 genome as BSgenome object
# path = "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/"
# setwd(dir = path)
# 
# forgeBSgenomeDataPkg(x = "BSgenome.seed", seqs_srcdir = ".")
## run following commands on CMD
# R.exe CMD build BSgenome.Afumigatus.AspGD.Af293
# R.exe CMD check BSgenome.Afumigatus.AspGD.Af293_03.05.06.tar.gz
# R.exe CMD INSTALL BSgenome.Afumigatus.AspGD.Af293_03.05.06.tar.gz










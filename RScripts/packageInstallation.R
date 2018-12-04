
rm(list = ls())

############################################################################
###############         uninstall packages          ########################
############################################################################

# create a list of all installed packages
ip <- as.data.frame(installed.packages())
head(ip)

# if you use MRO, make sure that no packages in this library will be removed
ip <- subset(ip, !grepl("MRO", ip$LibPath))

# we don't want to remove base or recommended packages either\
ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]

# determine the library where the packages are installed
path.lib <- unique(ip$LibPath)

# create a vector with all the names of the packages you want to remove
pkgs.to.remove <- ip[,1]
head(pkgs.to.remove)

# remove the packages
sapply(pkgs.to.remove, remove.packages, lib = path.lib)


############################################################################
###############         install packages            ########################
############################################################################

rm(list = ls())

install.packages(c("installr", "Rcpp", "tidyverse", "ggplot2", "ggpubr", "dplyr",
                   "tidyr", "data.table", "tibble", "purrr", "stringr", "readxl",
                   "lazyeval", "dendsort", "dendextend", "dynamicTreeCut", "RColorBrewer",
                   "hashmap", "reshape", "FactoMineR", "factoextra", "VennDiagram",
                   "imputeTS", "summarytools", "UpSetR", "esquisse", "corrgram"), dependencies = T)


if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("BiocGenerics", "IRanges", "GenomicRanges", "BiocGenerics", "S4Vectors", "rtracklayer",
           "Rsamtools", "BSgenome", "BiocParallel", "GenomeInfoDbData", "DESeq2", "ballgown",
           "pathview", "DO.db", "clusterProfiler", "topGO", "GenomicAlignments", "GenomicFeatures",
           "DiffBind", "ChIPComp", "regioneR", "tximport", "apeglm", "AnnotationHub", "Gviz",
           "KEGGREST", "KEGG.db", "KEGGprofile", "preprocessCore"))


install.packages(c("devtools"), dependencies = T)
library(devtools)
devtools::install_github("jokergoo/circlize")
# devtools::install_github("jokergoo/ComplexHeatmap")
# devtools::install_github("jokergoo/EnrichedHeatmap")

BiocManager::install("ComplexHeatmap")
BiocManager::install("EnrichedHeatmap")



#Time series data imputation
install.packages("imputeTS", dependencies = T)


## install my own packages
library(devtools)
devtools::document(pkg = "E:/Chris_UM/GitHub/chipmine")
devtools::install("E:/Chris_UM/GitHub/chipmine")
# devtools::install_github("lakhanp1/chipmine")

##########################################################################
###############         Org.Db packages            #######################
##########################################################################
## A. nidulans
library(devtools)
devtools::install_github("lakhanp1/fungal_resources/A_nidulans/org.Anidulans.eg.db")

## C. albicans
library(devtools)
devtools::install_github("lakhanp1/fungal_resources/C_albicans/org.Calbicans.eg.db")

## C. auris
library(devtools)
devtools::install_github("lakhanp1/fungal_resources/C_auris/org.Cauris.eg.db")


############################################################################
###############         BSgenome packages            #######################
############################################################################

## install A_fumigatus_293 genome as BSgenome object
path = "E:/Chris_UM/Database/A_fumigatus_293_version_s03-m05-r06/"
setwd(dir = path)

forgeBSgenomeDataPkg(x = "BSgenome.seed", seqs_srcdir = ".")
## run following commands on CMD
# R.exe CMD build BSgenome.Afumigatus.AspGD.Af293
# R.exe CMD check BSgenome.Afumigatus.AspGD.Af293_03.05.06.tar.gz
# R.exe CMD INSTALL BSgenome.Afumigatus.AspGD.Af293_03.05.06.tar.gz

library(BSgenome.Afumigatus.AspGD.Af293)











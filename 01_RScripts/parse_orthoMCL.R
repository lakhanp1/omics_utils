library(dplyr)
library(data.table)
library(tibble)
library(tidyr)



rm(list = ls())


path = "E:/Chris_UM/Database/Candida_auris_B8441/functionalAnnotation"
setwd(path)


dt = fread(file = "orthoMCL_cluster_report_ds8genomes", sep = "\t", header = T, na.strings = c(""), stringsAsFactors = F, select = c("fam", "Candida auris (B8441 - i)", "Candida albicans")) %>%
  dplyr::rename(CAuris = `Candida auris (B8441 - i)`, CAlbicans = `Candida albicans`)

df = dt %>% dplyr::filter(!is.na(CAuris) & !is.na(CAlbicans)) %>%
  tidyr::unnest(CAuris = strsplit(CAuris, " ")) %>%
  tidyr::unnest(CAlbicans = strsplit(CAlbicans, " "))







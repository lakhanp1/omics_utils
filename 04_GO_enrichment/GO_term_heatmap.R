
library(data.table)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

## This script plots a heatmap of gene's association with various GO SLIM terms


rm(list = ls())

path = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/HSP90_nComms_diff"
setwd(path)


gene2GoMatrixFile = "E:/Chris_UM/Database/GO_associations/goSlimMatrix.CAlbicans.tab"
geneFile = "Hsp90_DEGs.tab"
otherAnFile = "Hap43_targets.txt"


selectGo = c("GO:0016192|vesicle-mediated transport",
             "GO:0042221|response to chemical",
             "GO:0042254|ribosome biogenesis",
             "GO:0042493|response to drug",
             "GO:0050789|regulation of biological process",
             "GO:0006810|transport",
             "GO:0016070|RNA metabolic process",
             "GO:0006950|response to stress")



newColNames = gsub(pattern = "GO:\\d+\\|", replacement = "", x = selectGo, perl = T)

renameCols = structure(selectGo, names = newColNames)

geneToGoData = fread(file = gene2GoMatrixFile, sep = "\t", header = T, stringsAsFactors = F) %>%
  dplyr::select(geneId, !!selectGo)



genesDf = fread(file = geneFile, sep = "\t", header = T, stringsAsFactors = F)
otherData = fread(file = otherAnFile, sep = "\t", header = T, stringsAsFactors = F)


data = dplyr::left_join(x = genesDf, y = geneToGoData, by = c("geneId" = "geneId")) %>%
  dplyr::left_join(y = otherData, by = c("geneId" = "geneId")) %>%
  dplyr::mutate_at(.vars = vars(!!!selectGo), .funs = funs(ifelse(is.na(.), 0, .))) %>%
  tibble::column_to_rownames(var = "geneId") %>%
  dplyr::rename(!!!renameCols)
  


anDf = data[c("upOrDown", "Hap43_targets")]
anLables = structure(c("DEG type", "Hap43 target"), names = colnames(anDf))

htMat = data[, newColNames]


ht = Heatmap(matrix = htMat,
             column_title = "GO SLIM association of C. albicans Hsp90 depletion DEGs",
             row_title = "DEGs",
             col = structure(c("#ff7f00", "grey95"), names = c(1, 0)),
             na_col = "grey95",
             split = anDf["upOrDown"],
             show_row_names = FALSE,
             heatmap_legend_param = list(title = "GO SLIM association", 
                                         at = c(1, 0), 
                                         labels = c("Yes", "No")),
             row_title_side = "left"
)



an = HeatmapAnnotation(df = anDf, which = "row",
                       gap = unit(1, "mm"),
                       col = list(upOrDown = c("up" = "red", "down" = "blue"),
                                  Hap43_targets = c("1" = "#a65628")),
                       na_col = "grey95",
                       annotation_legend_param = list(
                         upOrDown = list(title = "DEG type"),
                         Hap43_targets = list(title = "Hap43 target",
                                              at = c(1), labels = c("Yes"))
                       )
)


htList = an + ht

png(filename = "Hsp90_DEGs_GO_SLIM_heatmap.png", width=3000, height=5000, res = 400)

draw(object = htList,
     row_dend_side = "left", heatmap_legend_side = "right", row_sub_title_side = "left")


## add column names
for(an in colnames(anDf)){
  decorate_annotation(annotation = an,
                      code = {
                        grid.text(label = anLables[[an]], x = unit(0.5, "npc"), y = unit(0, "npc") - unit(2, "mm"),
                                  default.units = "npc", just = "right", rot = 90)},
                      slice = 2
  )
}



dev.off()



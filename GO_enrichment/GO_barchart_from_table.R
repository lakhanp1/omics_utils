library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)


rm(list = ls())


path <- "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/GO_analysis"
setwd(path)


###########################################################################
goFile = "CAlb_Hsp90dep_nonOrtho_GOSLIM.txt"

goData = fread(file = goFile, sep = "\t", stringsAsFactors = F, header = T)

goData = dplyr::filter(goData, DEG_type != "allDEGs" & GO_term != "other" & GO_term != "biological_process")
  


goTable = data.table::dcast(data = goData, formula = GO_id + GO_term ~ DEG_type, value.var = "percent") %>%
  dplyr::filter_at(.vars = c("upDEGs", "downDEGs"), any_vars(. != 0))


filteredGoData = data.table::melt(data = goTable, id.vars = c("GO_id", "GO_term"),
                                  measure.vars = c("upDEGs", "downDEGs"),
                                  variable.name = "DEG_type",
                                  value.name = "percent")

cols = c("upDEGs" = "#E69F00", "downDEGs" = "#56B4E9")

title = "C. albicans Hsp90 depletion DEGs which do not have C. auris ortholog"

p1 = ggplot(data = filteredGoData, mapping = aes(x = GO_term, y = percent)) +
  geom_bar(mapping = aes(fill = DEG_type), stat="identity", position=position_dodge()) +
  scale_fill_manual(
    name = "DEGs",
    breaks = c("upDEGs", "downDEGs"),
    labels = c("Up", "Down"),
    values = cols) +
  ggtitle(str_wrap(title, 80)) +
  ylab("Percent of genes") +
  xlab("GO SLIM term") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold")) + 
  coord_flip()


png(filename = "CAlb_Hsp90dep_nonOrtho_GOSLIM.png", width=5000, height=6000, res = 500)
print(p1)
dev.off()

###########################################################################




###########################################################################
## comparison of two GO slim tables

goFile2 = "CAlb_HSP90depletion_GO_SLIM_BP.txt"

albicandGoData = fread(file = goFile2, sep = "\t", stringsAsFactors = F, header = T)

albicandGoData = dplyr::filter(albicandGoData, DEG_type != "allDEGs")

title2 = "GO SLIM comparison between C. albicans and C. auris HSP90 depletion response"

p2 = ggplot() +
  geom_bar(data = goData, 
           mapping = aes(x = GO_term, y = percent, fill = DEG_type),
           stat="identity", position=position_dodge()) +
  geom_bar(data = albicandGoData,
           mapping = aes(x = GO_term, y = - percent, fill = DEG_type),
           stat="identity", position=position_dodge()) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Percent of genes", 
                     breaks = seq(-40, 30, 20), 
                     labels = abs(seq(-40, 30, 20))) +
  scale_fill_manual(
    name = "DEGs",
    breaks = c("upDEGs", "downDEGs"),
    labels = c("Up", "Down"),
    values = cols) +
  ggtitle(str_wrap(title2, 80)) +
  ylab("Percent of genes") +
  xlab("GO SLIM term") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold")) + 
  coord_flip()


png(filename = "CAuris_CAlbicans_hsp90depletion_GO_SLIM_comparison.png", width=6000, height=6000, res = 500)
print(p2)
dev.off()


###########################################################################














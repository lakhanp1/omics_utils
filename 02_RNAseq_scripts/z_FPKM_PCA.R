
library(data.table)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(tibble)
library(ggrepel)

rm(list = ls())

path = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/correlationAnalysis"
setwd(path)


caOrthFile = "DEGs_CAlbicansOrthologs.tab"
hsf1ExpFile = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/HSF1_data/HSF1_OE_WT_DEL.tab"
sraExpFile = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/HSF1_data/SRA_FPKM_mat.tab"
cAurisExpFile = "E:/Chris_UM/Analysis/26_Cowen_CAuris_RNAseq/stringTie/FPKM_matrix.tab"
sampleInfoFile = "sampleInfo.tab"


##################################################################
## prepare the data
caOrth = fread(file = caOrthFile, sep = "\t", header = T, stringsAsFactors = F)
hsf1Exp = fread(file = hsf1ExpFile, sep = "\t", header = T, stringsAsFactors = F)
sraExp = fread(file = sraExpFile, sep = "\t", header = T, stringsAsFactors = F)
cAurisExp = fread(file = cAurisExpFile, sep = "\t", header = T, stringsAsFactors = F)
sampleDetails = fread(file = sampleInfoFile, sep = "\t", header = T, stringsAsFactors = F)


mergedData = left_join(x = caOrth, y = cAurisExp, by = c("geneId" = "geneId")) %>%
  left_join(y = hsf1Exp, by = c("ca_orf19_id" = "tracking_id")) %>%
  left_join(y = sraExp, by = c("ca_orf19_id" = "gene")) %>%
  # dplyr::select(colnames(caOrth), ends_with("_meanFPKM")) %>%
  mutate_if(.predicate = is.numeric, .funs = funs(log2(. + 0.001))) %>%
  dplyr::select(-ca_orf19_id)

# testData = head(mergedData)
# dcast.data.table(melt.data.table(testData, id.vars = "geneId"), variable ~ geneId)

exprData = dcast.data.table(melt.data.table(as.data.table(mergedData), id.vars = "geneId"), variable ~ geneId) %>%
  dplyr::rename(sampleId = variable) %>%
  dplyr::mutate(sampleId = as.character(sampleId)) %>%
  left_join(y = sampleDetails, by = c("sampleId" = "sampleId")) %>%
  dplyr::select(sampleId, project, library_name, condition, treatment, strain, everything()) %>%
  filter(!is.na(project))


projects = c("PRJNA167842", "Calbicans_HSF1")

subData = base::as.data.frame(exprData)  %>%
  filter(project %in% c(!!!projects)) %>%
  column_to_rownames(var = "sampleId")


##################################################################
## PCA
## tutorial followed: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#basics

res.pca <- PCA(subData, graph = FALSE, scale.unit = TRUE, quali.sup = 1:5)

eig.val <- get_eigenvalue(res.pca)

## scree plot: variance by PC
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))



##################################################################
## Graph of individuals
ind <- get_pca_ind(res.pca)

fviz_pca_ind(res.pca,
             # col.ind = exprData$treatment,
             fill.ind = subData$project,
             pointshape = 21,
             repel = TRUE,
             mean.point = FALSE,
             legend.title = "Study",
             pointsize = 3
)
  
## prepare the plot dataframe for ggplot
plotData = as.data.frame(ind$coord) %>%
  rownames_to_column(var = "sampleId") %>%
  left_join(y = sampleDetails, by = c("sampleId" = "sampleId"))

pltTitle = paste("PCA of samples from projects", paste(projects, collapse = " and "), sep = " ")
pointCol = structure(c("#E69F00", "#56B4E9"), names = projects)

pcaPlot = ggplot(data = plotData, mapping = aes(x = Dim.1, y = Dim.2)) +
  geom_point(mapping = aes(color = project, shape = treatment), size = 3) +
  geom_text_repel(mapping = aes(label = library_name), size = 3) +
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = pointCol) +
  xlab( paste("PC1 (", sprintf("%.2f", eig.val[1, "variance.percent"]), "%)", sep = "") ) +
  ylab( paste("PC2 (", sprintf("%.2f", eig.val[2, "variance.percent"]), "%)", sep = "") ) +
  ggtitle(pltTitle) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(face = "bold"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"))


##################################################################
##function to do PCA for each combination of experiment data


getPcaPlot = function(df, prj){
  projects = c("CAuris_HSP90", prj)
  
  subData = base::as.data.frame(df)  %>%
    filter(project %in% c(!!!projects)) %>%
    column_to_rownames(var = "sampleId")
  
  res.pca <- PCA(subData, graph = FALSE, scale.unit = TRUE, quali.sup = 1:5)
  
  eig.val <- get_eigenvalue(res.pca)
  
  ind <- get_pca_ind(res.pca)
  
  ## prepare the plot dataframe for ggplot
  plotData = as.data.frame(ind$coord) %>%
    rownames_to_column(var = "sampleId") %>%
    left_join(y = sampleDetails, by = c("sampleId" = "sampleId"))
  
  pltTitle = paste("PCA of samples from projects", paste(projects, collapse = " and "), sep = " ")
  pointCol = structure(c("#E69F00", "#56B4E9"), names = projects)
  
  pcaPlot = ggplot(data = plotData, mapping = aes(x = Dim.1, y = Dim.2)) +
    geom_point(mapping = aes(color = project, shape = treatment), size = 4) +
    geom_text_repel(mapping = aes(label = library_name), size = 4) +
    geom_hline(yintercept = 0, linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = pointCol) +
    xlab( paste("PC1 (", sprintf("%.2f", eig.val[1, "variance.percent"]), "%)", sep = "") ) +
    ylab( paste("PC2 (", sprintf("%.2f", eig.val[2, "variance.percent"]), "%)", sep = "") ) +
    ggtitle(pltTitle) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(face = "bold"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13, face = "bold"))
  
  png(filename = paste(prj, "_meanProfiles.png", sep = ""), width = 4000, height = 4000, res = 420)
  print(pcaPlot)
  dev.off()
  
  return(1)
  # return(pcaPlot)
}


plotList = list()
plt = getPcaPlot(df = exprData, prj = "PRJNA356057")


exprData %>%
  filter(project != "CAuris_HSP90") %>%
  group_by(project) %>%
  dplyr::do(pt = getPcaPlot(df = exprData, prj = unique(.$project)))











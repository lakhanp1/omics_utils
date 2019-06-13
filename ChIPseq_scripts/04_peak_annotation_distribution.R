library(chipmine)
library(here)

## plot peak annotation % distribution bar plot for multiple samples
rm(list = ls())

outDir <- here::here("analysis", "ChIP_summary")


if(!dir.exists(outDir)){
  dir.create(path = outDir)
}


file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")

TF_dataPath <- here::here("data", "ChIPseq_data")

outPrefix <- paste(outDir, "/CTCF_ChIP", sep = "")

##################################################################################
## get the sample details
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   dataPath = TF_dataPath)

exptDataList <- purrr::transpose(exptData) %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", exptData$sampleId, sep = ""), names = exptData$sampleId) },
  simplify = F, USE.NAMES = T
)

##################################################################################


peakTypeCounts <- tibble::tibble(sampleId = character(),
                                 peakType = character(),
                                 count = numeric())

i <- 1

for (i in 1:nrow(exptData)) {
  peakAn <- import_peak_annotation(sampleId = exptData$sampleId[i],
                                   peakAnnoFile = exptData$narrowpeakAnno[i])
  
  
  peakTypeTable <- table(peakAn[[tfCols$peakType[exptData$sampleId[i]]]])
  
  peakTypeCounts <- dplyr::bind_rows(
    peakTypeCounts,
    tibble(sampleId = exptData$sampleId[i],
           peakType = names(peakTypeTable),
           count = peakTypeTable,
           percent = peakTypeTable/sum(peakTypeTable))
  )
}

peakTypeCounts$peakType <- factor(
  peakTypeCounts$peakType,
  levels = c("upstream", "promoter", "5UTR", "tx_start", "include_tx", "EXON",
             "INTRON", "3UTR", "tx_end", "intergenic")
  )

peakTypeColors <- structure(
  RColorBrewer::brewer.pal(n = 10, name = "Paired"),
  names = c("upstream", "promoter", "5UTR", "tx_start", "include_tx", "EXON",
                     "INTRON", "3UTR", "tx_end", "intergenic"))

pt <- ggplot(data = peakTypeCounts) +
  geom_bar(mapping = aes(x = sampleId, y = count, fill = peakType),
           stat="identity", position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expand_scale(add = 0.01)) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(name = "CTCF peak position", values = peakTypeColors) +
  labs(title = "CTCF peak position distribution across genomic features") +
  theme(
    axis.text = element_text(size = 16),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16),
    title = element_text(size = 18),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    panel.grid = element_blank()
  )


# pdf(file = paste(outPrefix, ".peakType.distribution.pdf", sep = ""), width = 15, height = 10)
png(filename = paste(outPrefix, ".peakType.distribution.png", sep = ""), width = 4500, height = 3000, res = 350)
plot(pt)
dev.off()






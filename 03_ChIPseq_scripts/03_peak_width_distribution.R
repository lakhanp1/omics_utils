library(chipmine)
library(here)
library(ggridges)

## peak width density distribution plot for multiple samples
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

##################################################################################
i <- 1

peakWd <- tibble::tibble(sampleId = character(), peakWidth = numeric())

for (i in 1:nrow(exptData)) {
  peaks <- import_peaks_as_df(file = exptData$narrowpeakFile[i],
                              sampleId = exptData$sampleId[i],
                              peakFormat = "narrowPeak",
                              peakCols = c("peakChr", "peakStart", "peakEnd", "peakId",
                                           "peakEnrichment", "peakPval"),
                              rename = FALSE)
  
  peakWd <- dplyr::bind_rows(
    peakWd,
    tibble::tibble(
      sampleId = exptData$sampleId[i],
      condition = exptData$condition[i],
      peakWidth = peaks$peakEnd - peaks$peakStart
    )
  )
  
}

conditionColor <- structure(
  RColorBrewer::brewer.pal(n = length(unique(exptData$condition)), name = "Set1"),
  names = unique(exptData$condition))

countSummary <- dplyr::group_by(peakWd, sampleId) %>% 
  dplyr::summarise(
    peakCount = paste("n=", n(), sep = ""),
    min = min(peakWidth),
    max = max(peakWidth),
    # quantile.01 = quantile(peakWidth, 0.01),
    # quantile.10 = quantile(peakWidth, 0.1),
    # quantile.25 = quantile(peakWidth, 0.25),
    # quantile.50 = quantile(peakWidth, 0.50),
    # quantile.75 = quantile(peakWidth, 0.75),
    # quantile.90 = quantile(peakWidth, 0.90),
    # quantile.99 = quantile(peakWidth, 0.99),
    ) %>% 
  dplyr::mutate(
    summaryText = paste(peakCount, ", min=",min,", max=", max, sep = "")
  )


pt <- ggplot(data = peakWd,
       mapping = aes(x = peakWidth, y = sampleId)) +
  stat_density_ridges(
    mapping = aes(fill = condition),
    alpha = 0.7, size = 0.5,
    # rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles = c(0.25, 0.5, 0.75),
    vline_size = 1, vline_alpha = 1
  ) +
  # geom_hline(mapping = aes(yintercept = sampleId)) +
  geom_text(data = countSummary,
            mapping = aes(x = 800, y = sampleId, label = summaryText),
            hjust = 0, vjust = 0, nudge_y = 0.25) +
  scale_fill_manual(
    name = "Cell line",
    values = conditionColor
  ) +
  labs(title = "CTCF peak width distribution",
       subtitle = "25%, 50% and 75% quantiles marked with vertical lines\ndistribution outlier tails removed by limitting x-axis to 1500",
       x = "CTCF peak width") +
  scale_x_continuous(expand = c(0, 0.3), limits = c(100, 1500),
                     breaks = c(200, 500, 1000, 1400)) +
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1))) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.title.y = element_blank(),
    title = element_text(size = 18),
    axis.text.y = element_text(size = 14, vjust = 0),
    axis.text.x = element_text(size = 14),
    panel.grid = element_blank()
  )


# pdf(file = paste(outPrefix, ".peakWidth.pdf", sep = ""), width = 10, height = 10)
png(filename = paste(outPrefix, ".peakWidth.png", sep = ""), width = 3000, height = 3000, res = 350)
plot(pt)
dev.off()

##################################################################################







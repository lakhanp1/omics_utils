library(tidyverse)
library(ggbeeswarm)

setwd("E:/Chris_UM/Database/A_flavus_NRRL3557/A_flavus_NRRL3357.NCBI/repeatMasker")

df <- suppressMessages(readr::read_tsv(file = "GCA_009017415.1_ASM901741v1_repeatMasker.tab"))


df <- df %>% 
  dplyr::mutate(
    length = abs(query_end - query_start),
    type = dplyr::case_when(
      repeat_class == "Simple_repeat" ~ repeat_class,
      repeat_class == "Low_complexity" ~ repeat_class,
      repeat_class == "rRNA" ~ repeat_class,
      repeat_class == "LTR/Gypsy" ~ paste(repeat_class, "\n", matching_repeat)
    )
  )


pt <- ggplot(data = df, mapping = aes(x = type, y = length)) +
  ggbeeswarm::geom_quasirandom(size = 2) +
  labs(title = "RepeatMasker repeat length distribution in A. flavus NRRL3357 genome") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_blank(),
    plot.title = element_text(size = 16, face = "bold")
  )

png(filename = "repeat_length_distribution.png", width = 2000, height = 2000, res = 240)
pt
dev.off()


dplyr::group_by(df, repeat_class, type) %>% 
  dplyr::summarise(
    n = n(),
    mean_len = mean(length),
    length_sum = sum(length)
  )








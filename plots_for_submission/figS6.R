library(tidyverse)
library(scales)

pcs <- seq(0, 100, 10)

egenes <- 
    "../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/1-map_eqtls/th_50/2-permutations/results/permutations_%d.significant.txt" %>%
    sprintf(pcs) %>%
    setNames(pcs) %>%
    map_df(~read_delim(., delim = " ", col_names = FALSE), .id = "f") %>%
    count(f) %>%
    mutate(f = as.integer(f)) %>%
    arrange(f)

tiff("./plots/S6_fig.tiff", width = 5, height = 3, units = "in", res = 300)
ggplot(egenes, aes(f, n)) + 
    geom_point(size = 2.5) + 
    geom_line() +
    scale_x_continuous(breaks = egenes$f) +
    scale_y_continuous(labels = comma) +
    scale_color_manual(values = cols) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Times")) +
    labs(x = "Number of PCs", y = "Number of eGenes")
dev.off()
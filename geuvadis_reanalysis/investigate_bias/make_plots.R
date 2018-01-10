library(tidyverse)

df <- read_tsv("./integrated_data.tsv")

png("./plots/expression_by_refallele.png", width = 10, height = 5, units = "in", res = 200)
ggplot() +
    geom_jitter(data = df, aes(genotype, resid, color = dist),
                width = .25, alpha = 1/2, size = .75) +
    scale_color_gradient(low = "white", high = "red") +
    geom_text(data = distinct(df, index, gene_name, rsid),
              aes(x = 1.75, y = 3, label = rsid), size = 3, color = "white") +
    facet_grid(index~gene_name, scales = "free") +
    theme_dark() +
    labs(color = "divergence from reference")
dev.off()

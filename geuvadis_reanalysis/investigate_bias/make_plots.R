library(tidyverse)

df <- read_tsv("./integrated_data.tsv")

png("./plots/expression_by_refallele.png", width = 10, height = 5, units = "in", res = 200)
ggplot() +
    geom_jitter(data = df, aes(genotype, resid, color = dist),
                width = .25, size = .75) +
    scale_color_gradient(low = "cornflowerblue", high = "darkred") +
    geom_text(data = distinct(df, index, gene_name, rsid),
              aes(x = 1.75, y = 3, label = rsid), size = 3) +
    facet_grid(index~gene_name, scales = "free") +
    labs(y = "expression", color = "divergence from reference")
dev.off()

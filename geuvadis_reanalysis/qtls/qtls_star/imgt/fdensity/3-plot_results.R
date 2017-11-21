library(tidyverse)

fenrich_out <- 
    read_delim("./density_tf_around_qtl.txt", col_names = FALSE, delim = " ") %>%
    mutate(pos = (X1 + X2)/2) %>%
    select(pos, n = X3)

png("../../../plots/fdensity.png", width = 4, height = 4, units = "in", res = 200)
ggplot(fenrich_out, aes(pos, n)) +
    geom_line() +
    labs(x = "Distance to eQTLs", y = "Number of annotations / kb") +
    theme_bw()
dev.off()

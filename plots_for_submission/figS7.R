devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(scales)
library(cowplot)

make_pca_plot <- function(PC_x, PC_y) {
    
    ggplot(pcs, aes_string(PC_x, PC_y)) +
        geom_point(aes(color = pop)) +
        ggsci::scale_color_npg() +
        scale_x_continuous(breaks = pretty_breaks(3)) +
        scale_y_continuous(breaks = pretty_breaks(3)) +
        guides(color = guide_legend(override.aes = list(size = 2))) +
        theme(text = element_text(size = 11, family = "Arial")) +
        labs(color = "Population")
}

geuvadis_pops <- select(geuvadis_info, subject = name, pop)

pcs <- 
    read_delim("../geuvadis_reanalysis/data/pca_genotypes/eur.pca", delim = " ") %>%
    mutate(SampleID = sub("^.+_(PC\\d+)$", "\\1", SampleID)) %>%
    filter(SampleID %in% paste0("PC", 1:100)) %>%
    gather(subject, value, -1) %>%
    spread(SampleID, value) %>%
    mutate_at(vars(-subject), function(x) x/sqrt(sum(x^2))) %>%
    inner_join(geuvadis_pops, by = "subject")


p1 <- make_pca_plot("PC1", "PC2")
p2 <- make_pca_plot("PC2", "PC3")
p3 <- make_pca_plot("PC3", "PC4")
p4 <- make_pca_plot("PC4", "PC5")

leg <- get_legend(p1)


tiff("./plots/S7_fig.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_grid(p1 + guides(color=FALSE), 
          p2 + guides(color=FALSE), 
          leg,
          p3 + guides(color=FALSE), 
          p4 + guides(color=FALSE), 
          nrow = 2, ncol = 3, rel_widths = c(3, 3, 1))
dev.off()
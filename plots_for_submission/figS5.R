devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(ggsci)

hla_genes <- gencode_hla %>%
    filter(! gene_name %in% c("HLA-DRA", "DPA1"))

caveman <- 
    "../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/caveman/results.hla" %>%
    read_tsv() %>%
    mutate(gene = factor(gene, levels = hla_genes$gene_name),
           rank = as.character(rank))

tiff("./plots/S5_fig.tiff", width = 7.5, height = 3, units = "in", res = 300)
ggplot(data = caveman, aes(rank, Probability, fill = index)) +
    geom_bar(stat = "identity", position = "dodge", alpha = .8) +
    scale_fill_npg() +
    facet_wrap(~gene, nrow = 1, scales = "free_x") +
    labs(x = "") +
    theme(text = element_text(size = 11, family = "Arial"), 
          legend.position = "top")
dev.off()

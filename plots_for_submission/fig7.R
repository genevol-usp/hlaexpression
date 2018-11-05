library(tidyverse)
devtools::load_all("/home/vitor/Libraries/hlaseqlib")

calc_ase <- function(tpm) min(tpm)/sum(tpm)

ase_df <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    group_by(subject, locus) %>%
    filter(n_distinct(allele) == 2) %>%
    summarise(ase = calc_ase(tpm)) %>%
    ungroup()

tiff("./plots/Fig7.tiff", width = 6, height = 3, units = "in", res = 300)
ggplot(ase_df, aes(reorder(locus, ase, median), ase)) +
    ggbeeswarm::geom_quasirandom(varwidth = TRUE, size = 1, alpha = 1/2) +
    stat_summary(fun.y = median, geom = "point", shape = "\U2014", size = 20, color = "#DC0000B2") +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(x = "", y = "ASE") +
    theme_bw() +
    theme(text = element_text(size = 10, family = "Arial"))
dev.off()
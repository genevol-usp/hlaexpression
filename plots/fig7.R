library(tidyverse)
devtools::load_all("/home/vitor/Libraries/hlaseqlib")

calc_ase <- function(tpm) min(tpm)/sum(tpm)

hla_loci <- paste0("HLA-", c("A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1", "DRA"))

ase_df <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_loci) %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    group_by(subject, locus) %>%
    filter(n_distinct(allele) == 2) %>%
    summarise(ase = calc_ase(tpm)) %>%
    ungroup()

tiff("./plots/Fig7.tiff", width = 12, height = 5, units = "cm", res = 300)
ggplot(ase_df, aes(reorder(locus, ase, median), ase)) +
    ggbeeswarm::geom_quasirandom(varwidth = TRUE, size = 1, alpha = 1/2) +
    stat_summary(fun.y = median, geom = "point", shape = "\U2014", size = 10, color = "#DC0000B2") +
    coord_cartesian(ylim = c(0, 0.5)) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Times")) +
    labs(x = NULL, y = "ASE")
dev.off()

ase_df %>% 
    group_by(locus) %>%
    summarise(e = mean(ase >= 0.4)) %>%
    ungroup() %>%
    summarise(round(mean(e) * 100))

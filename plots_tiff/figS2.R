devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(cowplot)
library(scales)


kallisto <- 
    "../geuvadis_reanalysis/expression/5-pseudoalignment/hla_personalized/quantifications_2/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup() %>%
    gather(unit, estimate, est_counts, tpm)

star_salmon <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup() %>%
    gather(unit, estimate, est_counts, tpm)

star_kallisto_df <- 
    left_join(star_salmon, kallisto, by = c("subject", "locus", "unit"),
              suffix = c(".mapping", ".pseudo")) %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name))

cor_df_star_kallisto <- star_kallisto_df %>%
    group_by(locus, unit) %>%
    do(data.frame(r = cor(.$estimate.mapping, .$estimate.pseudo),
                  p = cor(.$estimate.mapping, .$estimate.pseudo, method = "spearman"),
                  x = min(.$estimate.mapping),
                  y = max(.$estimate.pseudo))) %>%
    ungroup() %>%
    mutate_at(vars(r, p), ~round(., digits = 2)) %>%
    mutate(label = paste("r == ", r, "*','~rho ==", p))

p_counts <- 
    ggplot(filter(star_kallisto_df, unit == "est_counts"),
           aes(estimate.mapping, estimate.pseudo)) +
    geom_abline() +
    geom_point(size = .8) +
    scale_x_continuous(breaks = pretty_breaks(n = 2), labels = comma) +
    scale_y_continuous(breaks = pretty_breaks(n = 3), labels = comma) +
    facet_wrap(~locus, scales = "free") +
    geom_text(data = filter(cor_df_star_kallisto, unit == "est_counts"), 
              aes(x, y, label = label), family = "Times",
              parse = TRUE, hjust = "inward", vjust = "inward", size = 3.5) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Times"),
          axis.text = element_text(hjust = 1)) +
    labs(x = "STAR-Salmon", y = "kallisto", title = "Estimated Counts")

p_tpm <- 
    ggplot(filter(star_kallisto_df, unit == "tpm"),
           aes(estimate.mapping, estimate.pseudo)) +
    geom_abline() +
    geom_point(size = .8) +
    scale_x_continuous(breaks = pretty_breaks(n = 2), labels = comma) +
    scale_y_continuous(breaks = pretty_breaks(n = 3), labels = comma) +
    facet_wrap(~locus, scales = "free") +
    geom_text(data = filter(cor_df_star_kallisto, unit == "tpm"), 
              aes(x, y, label = label), family = "Times",
              parse = TRUE, hjust = "inward", vjust = "inward", size = 3.5) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Times"),
          axis.text = element_text(hjust = 1)) +
    labs(x = "STAR-Salmon", y = "kallisto", title = "Transcripts per Million")



tiff("./plots/S2_fig.tiff", width = 6, height = 8, units = "in", res = 300)
plot_grid(p_counts, NULL, p_tpm, nrow = 3, ncol = 1, rel_heights = c(1, 0.07, 1))
dev.off()

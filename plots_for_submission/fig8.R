devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(cowplot)
library(GGally)
#library(scales)
#library(ggpmisc)
#library(ggrepel)

calc_trans_cors <- function(locus1, locus2, df) {
    
    m_sub <- df %>% 
        select(subject, locus1, locus2) %>%
        group_by(subject) %>%
        mutate_at(vars(locus2), rev) %>%
        ungroup() %>%
        select(-subject) %>%
        cor(use = "pairwise.complete.obs")
    
    m[c(locus1, locus2), c(locus1, locus2)] <<- m_sub 
}

hla_genes <- paste0("HLA-", c("A", "B", "C", "DPB1", "DRB1", "DQA1", "DQB1"))

hlapers <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    mutate(subject = convert_ena_ids(subject),
           allele = gsub("IMGT_", "", allele)) %>%
    select(subject, locus, allele, tpm)

hlapers_gene <- hlapers %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

global_cors <- hlapers_gene %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    spread(locus, tpm) %>% 
    select(-subject) %>% 
    ggcorr(label = TRUE, label_round = 2, label_size = 3.25, hjust = 0.6) +
    labs(title = "Gene-level") +
    theme(legend.position = "none",
          text = element_text(size = 9, family = "Arial"),
          plot.title = element_text(size = 11, family = "Arial", hjust = 0.5))


haps_data <- "../geuvadis_reanalysis/phase_hla/phase_hla_haps_snps.tsv" %>%
    read_tsv() %>% 
    distinct(subject, locus, hap, allele_gene) %>%
    rename(allele = allele_gene) %>%
    left_join(hlapers, by = c("subject", "locus", "allele")) %>%
    distinct() %>%
    select(-allele) %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    spread(locus, tpm)

cis_cors <- haps_data %>% 
    select(-subject, -hap) %>%
    ggcorr(label = TRUE, label_round = 2, label_size = 3.25, hjust = 0.6) +
    labs(title = "Within haplotypes") +
    theme(legend.position = "none",
          text = element_text(size = 9, family = "Arial"),
          plot.title = element_text(size = 11, family = "Arial", hjust = 0.5))

m <- 
    matrix(NA, nrow = 7, ncol = 7, 
           dimnames = list(sort(sub("HLA-", "", hla_genes)), 
                           sort(sub("HLA-", "", hla_genes))))

phase_data <- tibble(locus1 = hla_genes, locus2 = hla_genes) %>% 
    mutate_all(~sub("HLA-", "", .)) %>% 
    expand(locus1, locus2) %>% 
    filter(locus1 != locus2) %>%
    pmap(calc_trans_cors, haps_data) 

trans_cors <- 
    ggcorr(data = NULL, cor_matrix = m, label = TRUE, label_round = 2,  
           label_size = 3.2, hjust = 0.6) +
    scale_fill_gradient2(name = "r", limits = c(0, 1), breaks = c(0, 0.5, 1), 
                         low = "#3B9AB2", mid = "#EEEEEE", high = "#F21A00") +
    labs(title = "Between haplotypes") +
    theme(text = element_text(size = 9, family = "Arial"),
          plot.title = element_text(size = 11, family = "Arial", hjust = 0.5))

legend <- get_legend(trans_cors)

trans_cors <- trans_cors + theme(legend.position = "none")

plot_correlations <- 
    plot_grid(global_cors, cis_cors, trans_cors, legend, nrow = 1, 
              rel_widths = c(1, 1, 1, .25))


classII_and_CIITA <- gencode_chr_gene %>%
    filter(gene_name %in% c("HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "CIITA"))

class_2_trans_df <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/gene_quantifications.tsv" %>%
    read_tsv() %>%
    inner_join(classII_and_CIITA, by = "gene_id") %>%
    select(subject, locus = gene_name, tpm) %>%
    spread(locus, tpm) %>%
    gather(locus, tpm, -subject, -CIITA) %>%
    select(subject, locus, tpm, CIITA) %>%
    arrange(subject, locus)

cor_df <- class_2_trans_df %>%
    group_by(locus) %>%
    do(data.frame(r = cor(.$tpm, .$CIITA),
                  x = min(.$tpm),
                  y = max(.$CIITA))) %>%
    ungroup() %>%
    mutate(r = round(r, digits = 2))

plot_ciita <- ggplot(class_2_trans_df, aes(tpm, CIITA)) +
    geom_point(size = 1) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(2)) +
    geom_text(data = cor_df, aes(x, y, label = paste("r =", r)),
              hjust = "inward", vjust = "inward", size = 4) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Arial"),
          axis.text = element_text(size = 11, family = "Arial", hjust = 1),
          strip.text = element_text(size = 11, family = "Arial")) + 
    facet_wrap(~locus, nrow = 1, scales = "free") +
    labs(x = NULL)


tiff("./plots/Fig8.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_grid(plot_correlations, plot_ciita, nrow = 2, ncol = 1, 
          rel_heights = c(1, .6), labels = c("A", "B"), label_size = 12, 
          hjust = 0)
dev.off()

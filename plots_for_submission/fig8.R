devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(scales)
library(cowplot)
library(GGally)

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

hla_genes <- gencode_hla$gene_name

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
    select(!!! syms(sub("HLA-", "", gencode_hla$gene_name))) %>% 
    ggcorr(label = TRUE, label_round = 2, label_size = 1.8, hjust = 0.8, size = 2.5, family = "Times") +
    labs(title = "Gene-level") +
    theme(legend.position = "none",
          text = element_text(size = 5, family = "Times"),
          plot.title = element_text(size = 10, family = "Times", hjust = 0.5))


phased <- "../geuvadis_reanalysis/phase_hla/concordant_set.tsv" %>%
    read_tsv() %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    select(subject, locus, hap, tpm) %>%
    spread(locus, tpm) %>%
    select(subject, hap, !!! syms(sub("HLA-", "", hla_genes)))


cis_cors <- phased %>% 
    select(-subject, -hap) %>%
    ggcorr(label = TRUE, label_round = 2, label_size = 1.8, hjust = 0.8, size = 2.5, family = "Times") +
    labs(title = "Within haplotypes") +
    theme(legend.position = "none",
          text = element_text(size = 9, family = "Times"),
          plot.title = element_text(size = 9, family = "Times", hjust = 0.5))

m <- 
    matrix(NA, nrow = length(hla_genes), ncol = length(hla_genes), 
           dimnames = list(sub("HLA-", "", hla_genes), 
                           sub("HLA-", "", hla_genes)))

phase_data <- tibble(locus1 = hla_genes, locus2 = hla_genes) %>% 
    mutate_all(~sub("HLA-", "", .)) %>% 
    expand(locus1, locus2) %>% 
    filter(locus1 != locus2) %>%
    pmap(calc_trans_cors, phased) 

trans_cors <- 
    ggcorr(data = NULL, cor_matrix = m, label = TRUE, label_round = 2,  
           label_size = 1.8, hjust = 0.8, size = 2.5, family = "Times") +
    scale_fill_gradient2(name = "r", limits = c(0, 1), breaks = c(0, 0.5, 1), 
                         low = "#3B9AB2", mid = "#EEEEEE", high = "#F21A00") +
    labs(title = "Between haplotypes") +
    theme(text = element_text(size = 9, family = "Times"),
          plot.title = element_text(size = 9, family = "Times", hjust = 0.5))

legend <- get_legend(trans_cors)

trans_cors <- trans_cors + theme(legend.position = "none")

plot_correlations <- 
    plot_grid(global_cors, cis_cors, trans_cors, legend, NULL, nrow = 1, 
              rel_widths = c(1, 1, 1, .25, 0.05))

classII_genes <-  c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1")

classII_and_CIITA <- gencode_chr_gene %>%
    filter(gene_name %in% c(classII_genes, "CIITA")) %>%
    arrange(as.integer(chr), start) %>%
    mutate(gene_name = factor(gene_name, levels = .$gene_name))

class_2_trans_df <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/gene_quantifications.tsv" %>%
    read_tsv() %>%
    inner_join(classII_and_CIITA, by = "gene_id") %>%
    select(subject, locus = gene_name, tpm) %>%
    spread(locus, tpm) %>%
    gather(locus, tpm, -subject, -CIITA) %>%
    mutate(locus = factor(locus, levels = classII_and_CIITA$gene_name)) %>%
    select(subject, locus, tpm, CIITA)

cor_df <- class_2_trans_df %>%
    group_by(locus) %>%
    do(data.frame(r = cor(.$tpm, .$CIITA),
                  x = min(.$tpm),
                  y = max(.$CIITA))) %>%
    ungroup() %>%
    mutate(r = round(r, digits = 2))

plot_ciita <- ggplot(class_2_trans_df, aes(tpm, CIITA)) +
    geom_point(size = .5, alpha = .5) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = pretty_breaks(2), labels = comma) +
    geom_label(data = cor_df, aes(x, y, label = paste("r =", r)),
              hjust = "inward", vjust = "inward", size = 3, family = "Times",
              label.padding = unit(0.05, "lines"), label.size = NA, alpha = 0.4) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times"),
          axis.text.x = element_text(angle = 15, hjust = 0.8)) +
    facet_wrap(~locus, nrow = 1, scales = "free") +
    labs(x = NULL)


tiff("./plots/Fig8.tiff", width = 13.2, height = 9, units = "cm", res = 300)
plot_grid(plot_correlations, plot_ciita, nrow = 2, ncol = 1, 
          rel_heights = c(1, .5), labels = c("A", "B"), label_size = 9, 
          label_fontfamily = "Times",
          hjust = -1)
dev.off()

devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(scales)
library(cowplot)
library(GGally)


## CRD
gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene_name, start, end, strand)

crd <- "../LCL_ALL.chr6.subset.txt.gz" %>%
    read_delim(col_names = FALSE, delim = " ") %>%
    select(-X2, -X4, -X6, -X8, -X9, -X11) %>%
    filter(between(X3, min(gencode_hla_v19$start) - 1e6, max(gencode_hla_v19$end) + 1e6),
           between(X7, min(gencode_hla_v19$start) - 1e6, max(gencode_hla_v19$end) + 5e5)) %>%
    mutate(index = (X1 + X5)/2L,
           pos = (X3 + X7)/2L,
           y = X5 - X1,
           r2 = X10^2) %>%
    group_by(index) %>%
    mutate(pos_label = min(pos)) %>%
    ungroup()
    
gene_pos <- gencode_hla_v19 %>%
    mutate(gene_name = sub("HLA-", "", gene_name),
           pos = ifelse(strand == "+", start, end),
           closest = crd$index[map_dbl(pos, ~which.min(abs(. - crd$pos_label)))])

crd_plot <- ggplot() +
    geom_point(data = crd, 
               aes(x = index, y = y, alpha = r2), 
               color = "blue", size = .5) +
    scale_alpha_continuous(range = c(0, 1)) +
    scale_x_continuous(breaks = crd$index[seq(1, nrow(crd), 4e4)],
                       labels = round(crd$pos_label[seq(1, nrow(crd), 4e4)]/1e6, 1)) +
    geom_hline(yintercept = -15, size = 2, color = "grey", alpha = 1/2) +
    geom_segment(data = filter(gene_pos, strand == "+"),
                 aes(x = closest, xend = closest+7.5, y = -15, yend = -15),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed", ends = "last"),
                 alpha = 1/2) +
    geom_segment(data = filter(gene_pos, strand == "-"),
                 aes(x = closest, xend = closest-7.5, y = -15, yend = -15),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed", ends = "last"),
                 color = "red", alpha = 1/2) +
    geom_text(data = filter(gene_pos, strand == "+"),
              aes(x = closest, y = -5, label = gene_name),
              size = 3, hjust = 0) +
    ggrepel::geom_text_repel(data = filter(gene_pos, strand == "-"),
                             aes(x = closest, y = -25, label = gene_name),
                             size = 3, color = "red") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = c(.95, .8)) +
    labs(x = "Genomic position (Mb)", alpha = expression(~r^2)) +
    guides(alpha = guide_legend(override.aes = list(size = 2.5)))


### Global correlations
 
star_imgt_tpm <- 
    read_tsv("./expression/star/imgt/quantifications_2/processed_imgt_quants.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele, tpm) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

global_cors <- star_imgt_tpm %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    spread(locus, tpm) %>%
    select(-subject) %>% ggcorr(label = TRUE)


### within vs between haplotypes

make_data <- function(locus1, locus2, df) {
    
    cis <- select(df, subject, locus1, locus2) %>%
        drop_na() %>%
        rename(gene1 = !!locus1, gene2 = !!locus2)
    
    trans <- cis %>% 
        group_by(subject) %>%
        mutate_at(vars(gene2), rev) %>%
        ungroup()
    
    bind_rows(list("Within haplotypes" = cis, 
                   "Between haplotypes" = trans), .id = "level") %>%
        mutate(pair = paste(locus1, "vs", locus2),
               level = factor(level, levels = c("Within haplotypes", "Between haplotypes"))) %>%
        select(pair, everything())
}

plotphase <- function(phase_df) {
    
    phase_cor_df <- phase_df %>%
        group_by(pair, level) %>%
        summarize(r = cor(gene1, gene2),
                  x = min(gene1),
                  y = max(gene2)) %>%
        ungroup() %>%
        mutate(r = round(r, digits = 2))
    
    ggplot(phase_df, aes(gene1, gene2)) +
        geom_point(size = 1) +
        geom_smooth(method = lm, se = FALSE) + 
        scale_x_continuous(breaks = pretty_breaks(2)) +
        scale_y_continuous(breaks = pretty_breaks(2)) +
        geom_text(data = phase_cor_df,
                  aes(x, y, label = paste("r =", r)),
                  hjust = "inward", vjust = "inward", size = 4) +
        facet_grid(pair~level, scales = "free") +
        theme_bw() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_blank(),
              strip.text = element_text(size = 8))
}

concordant_haps_I <- 
    read_tsv("./expression/star/phase_hla_alleles/data/concordant_haps_classI.tsv") %>%
    select(subject, hap)

concordant_haps_II <- 
    read_tsv("./expression/star/phase_hla_alleles/data/concordant_haps_classII.tsv") %>%
    select(subject, hap)

tpm_by_allele_I <- 
    read_tsv("./expression/star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv") %>%
    select(subject, locus, hap, tpm) %>%
    filter(locus %in% c("A", "B", "C")) %>%
    spread(locus, tpm) %>%
    inner_join(concordant_haps_I, by = c("subject", "hap")) %>%
    arrange(subject, hap)

tpm_by_allele_II <- 
    read_tsv("./expression/star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv") %>%
    select(subject, locus, hap, tpm) %>%
    filter(!locus %in% c("A", "B", "C")) %>%
    spread(locus, tpm) %>%
    inner_join(concordant_haps_II, by = c("subject", "hap")) %>%
    arrange(subject, hap)

tpm_by_allele <- 
    full_join(tpm_by_allele_I, tpm_by_allele_II, by = c("subject", "hap"))

phase_data <- 
    tribble(
        ~locus1, ~locus2,
        "A"    , "B",
        "A"    , "C",
        "B"    , "C",
        "DQA1" , "DQB1",
        "DQA1" , "DRB1") %>%
    pmap_df(make_data, tpm_by_allele)

phase_list <- phase_data %>% split(.$pair)

p1 <- plotphase(phase_list[[1]])
p2 <- plotphase(phase_list[[2]]) + theme(strip.text.x = element_blank())    
p3 <- plotphase(phase_list[[3]]) + theme(strip.text.x = element_blank())   
p4 <- plotphase(phase_list[[4]]) + theme(strip.text.x = element_blank())   
p5 <- plotphase(phase_list[[5]]) + theme(strip.text.x = element_blank())

within_vs_between <- 
    plot_grid(p1, p2, p3, p4, p5, ncol = 1, rel_heights = c(1, .9, .9, .9, .9), labels = "C")

plot_AB <- plot_grid(crd_plot, global_cors, ncol = 1, labels = c("A", "B"), 
                     rel_heights = c(.55, .45))

png("./correlations.png", width = 10, height = 6.5, units = "in", res = 300)
plot_grid(plot_AB, within_vs_between, ncol = 2, rel_widths = c(.55, .45))
dev.off()

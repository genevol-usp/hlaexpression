devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(cowplot)
library(scales)


dist_to_ref <- "../imgt_index/distance_to_ref/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)

hlapers <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    mutate(subject = convert_ena_ids(subject),
           allele = gsub("IMGT_", "", allele)) %>%
    select(subject, locus, allele, tpm)


haps_expression <- "../geuvadis_reanalysis/phase_hla/phase/phase_hla_haps_snps.tsv" %>%
    read_tsv() %>%
    filter(rank == 0) %>% 
    group_by(subject, locus) %>%
    filter(!any(uncertain_gene == 1)) %>%
    ungroup() %>%
    left_join(hlapers, by = c("subject", "locus", "allele_gene" = "allele")) %>%
    distinct(subject, locus, hap, .keep_all = TRUE) %>%
    mutate(locus = sub("HLA-", "", locus),
           locus = factor(locus, levels = sub("HLA-", "", gencode_hla$gene_name)),
           lineage = hla_trimnames(allele_gene, 1)) %>%
    select(subject, hap, locus, allele = allele_gene, lineage, tpm, qtl_rsid = rsid, 
           qtl_allele = allele_snp)

hap_hla_genot <- haps_expression %>%
    select(subject, locus, hap, allele)

hap_snps <- haps_expression %>%
    select(subject, locus, hap, allele = qtl_allele)

phen_best <- 
    "../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/1-map_eqtls/th_50/1-phenotypes/phenotypes_60.bed.gz" %>%
    read_tsv() %>%
    inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
    select(gene_name, HG00096:NA20828) %>%
    gather(subject, resid, -gene_name) %>%
    select(subject, locus = gene_name, resid) %>%
    mutate(locus = sub("HLA-", "", locus),
           locus = factor(locus, levels = sub("HLA-", "", gencode_hla$gene_name)))

qtls_high_low <- left_join(hap_snps, phen_best, by = c("subject", "locus")) %>%
    group_by(locus, allele) %>% 
    summarize(eQTL = mean(resid)) %>%
    group_by(locus) %>%
    mutate(eQTL = ifelse(eQTL == max(eQTL), "High", "Low")) %>%
    ungroup()

lineage_phased <- haps_expression %>%
    left_join(qtls_high_low, by = c("locus", "qtl_allele" = "allele")) %>%
    select(subject, locus, allele, lineage, qtl_rsid, eQTL, tpm) %>%
    left_join(dist_to_ref, by = "allele") %>%
    mutate(eQTL = factor(eQTL, levels = c("Low", "High"))) %>%
    group_by(lineage) %>%
    filter(n() >= 10) %>%
    ungroup()

best_eqtl_locus <- distinct(lineage_phased, locus, variant = qtl_rsid)

eqtl_info <- "../geuvadis_reanalysis/phase_hla/phase/eqtl_snps.vcf" %>%
    read_tsv(comment = "##") %>%
    select(-`#CHROM`, -POS, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    gather(subject, genotype, -(ID:ALT)) %>%
    inner_join(best_eqtl_locus, by = c("ID" = "variant")) %>%
    separate(genotype, c("h1", "h2"), convert = TRUE) %>%
    gather(hap, allele, h1:h2) %>%
    mutate(allele = ifelse(allele == 0, REF, ALT)) %>%
    arrange(subject, locus, ID, allele) %>%
    group_by(subject, locus, ID) %>%
    summarize(genotype = paste(sort(allele), collapse = "/")) %>%
    select(subject, locus, variant = ID, genotype) %>%
    ungroup() %>%
    unite(id, genotype, locus, sep = "_", remove = FALSE)

eqtls_expression_df <- left_join(eqtl_info, phen_best, by = c("subject", "locus"))



plot_lineages <- ggplot(data = lineage_phased,
       aes(x = reorder(lineage, tpm, FUN = median, na.rm = TRUE), y = tpm)) +
    geom_jitter(aes(color = eQTL), size = .7) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "grey5", alpha = .1) +
    scale_color_manual(values = c("Low" = ggsci::pal_npg()(6)[6], 
                                  "High" = ggsci::pal_npg()(1), 
                                  "ND" = "grey")) +
    scale_x_discrete(labels = function(x) sub("^(.+\\*)", "", x)) +
    scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) +
    facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times"), 
          strip.text = element_text(face = "bold"), 
          legend.position = "top") +
    labs(x = "HLA lineage", y = "TPM")

legend <- get_legend(plot_lineages)

plot_lineages_final <- plot_lineages + theme(legend.position = "none")


plot_slopes <- ggplot(data = eqtls_expression_df, 
       aes(reorder(id, resid, "mean"), resid)) +
    geom_jitter(width = .25, alpha = 1/2, size = .75) +
    geom_smooth(aes(group = 1), method = lm, se = FALSE) +
    scale_x_discrete(labels = function(x) sub("^([^_]+).+$", "\\1", x)) +
    scale_y_continuous(position = "right") +
    geom_label(data = distinct(eqtls_expression_df, locus, variant), 
              aes(x = 0.5, y = 2.7, label = variant), 
              label.padding = unit(0.05, "lines"), label.size = NA, alpha = 0.4,
              size = 2.5, family = "Times", hjust = "inward") +
    coord_cartesian(ylim = c(-3, 3.2)) +
    facet_wrap(~locus, ncol = 1, scales = "free") +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Times"), 
          strip.text = element_blank()) +
    labs(x = "eQTL genotype", y = "PCA-corrected expression")

plot_grid <- plot_grid(legend, NULL, plot_lineages_final, plot_slopes, ncol = 2, 
                   rel_widths = c(2.5, 1), rel_heights = c(.07, 1))


tiff("./plots/Fig6.tiff", width = 11, height = 18, units = "cm", res = 300)
ggdraw(plot_grid)
dev.off()



lineage_df %>%
    select(subject, locus, lineage, tpm) %>%
    group_by(lineage) %>%
    filter(n() > 10) %>%
    group_by(locus) %>%
    filter(n_distinct(lineage) > 1) %>% 
    ungroup() %>% 
    split(.$locus) %>%
    map(~oneway.test(tpm~lineage, data = .) %>% broom::tidy()) %>%
    bind_rows(.id = "locus") %>%
    select(locus, `num df`, `denom df`, F = statistic, p.value) %>%
    write_tsv("./f_onewaytest_lineages.tsv")

lineage_df %>%
    select(subject, locus, lineage, tpm) %>%
    group_by(lineage) %>%
    filter(n() > 10) %>%
    group_by(locus) %>%
    filter(n_distinct(lineage) > 1) %>% 
    ungroup() %>% 
    split(.$locus) %>%
    map(~lm(tpm~lineage, data = .) %>% anova() %>% broom::tidy()) %>%
    bind_rows(.id = "locus") %>%
    filter(term == "lineage") %>%
    select(locus, df, F = statistic, p.value) %>%
    write_tsv("./f_test_lineages.tsv")



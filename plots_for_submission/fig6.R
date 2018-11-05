devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(cowplot)
#library(GGally)
library(scales)
#library(ggpmisc)
#library(ggrepel)

lineage_df <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele, tpm) %>%
    mutate(subject = convert_ena_ids(subject),
           allele = gsub("IMGT_", "", allele),
           allele_3f = hla_trimnames(allele, 3),
           allele_2f = hla_trimnames(allele, 2),
           lineage = hla_trimnames(allele, 1))

dist_to_ref <- "../imgt_index/distance_to_ref/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)

haps_expression <- "../geuvadis_reanalysis/phase_hla/phase_hla_haps_snps.tsv" %>%
    read_tsv() %>%
    filter(rank == 0L) %>%
    rename(hla_allele = allele_gene) %>%
    left_join(lineage_df, by = c("subject", "locus", "hla_allele" = "allele")) %>%
    distinct() 

hap_hla_genot <- haps_expression %>%
    select(subject, locus, hap, hla_allele)

hap_snps <- haps_expression %>%
    select(subject, locus, hap, allele = allele_snp)

phen_best <- 
    "../geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/1-map_eqtls/th_50/1-phenotypes/phenotypes_60.bed.gz" %>%
    read_tsv() %>%
    inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
    select(gene_name, HG00096:NA20828) %>%
    gather(subject, resid, -gene_name) %>%
    select(subject, locus = gene_name, resid)

qtls_high_low <- left_join(hap_snps, phen_best, by = c("subject", "locus")) %>%
    group_by(locus, allele) %>% 
    summarize(eQTL = mean(resid)) %>%
    group_by(locus) %>%
    mutate(eQTL = ifelse(eQTL == max(eQTL), "High", "Low")) %>%
    ungroup()

lineage_phased <- haps_expression %>%
    left_join(qtls_high_low, by = c("locus", "allele_snp" = "allele")) %>%
    select(subject, locus, hla_allele, lineage, eQTL, tpm) %>%
    left_join(dist_to_ref, by = c("hla_allele" = "allele")) %>%
    mutate(eQTL = factor(eQTL, levels = c("Low", "High"))) %>%
    group_by(lineage) %>%
    filter(n() >= 10L) %>%
    ungroup() %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name))

best_eqtl_locus <- read_tsv("../geuvadis_reanalysis/eqtl_mapping/plots/eqtl.tsv") %>%
    filter(rank == 0) %>%
    select(locus = gene, variant = var_id)

eqtl_info <- "../geuvadis_reanalysis/phase_hla/eqtl_snps.vcf" %>%
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

eqtls_expression_df <- left_join(eqtl_info, phen_best, by = c("subject", "locus")) %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name))


p1 <- ggplot(data = lineage_phased,
             aes(x = reorder(lineage, tpm, FUN = median, na.rm = TRUE), 
                 y = tpm)) +
    geom_jitter(aes(color = eQTL), size = .75) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "grey5", alpha = .1) +
    scale_color_manual(values = c("Low" = ggsci::pal_npg()(6)[6], 
                                  "High" = ggsci::pal_npg()(1), 
                                  "ND" = "grey")) +
    scale_x_discrete(labels = function(x) sub("^(.+\\*)", "", x)) +
    scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) +
    facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
    labs(x = " ", y = "TPM") +
    theme_bw() +
    theme(text = element_text(size = 10, family = "Arial"),
          strip.text = element_text(face = "bold"),
          legend.position = "top") +
    guides(color = guide_legend(override.aes = list(size = 4)))

legend <- get_legend(p1)

p1 <- p1 + theme(legend.position = "none")

p2 <- ggplot(eqtls_expression_df, aes(reorder(id, resid, "mean"), resid)) +
    geom_jitter(width = .25, alpha = 1/2, size = .75) +
    geom_smooth(aes(group = 1), method = lm, se = FALSE) +
    scale_x_discrete(labels = function(x) sub("^([^_]+).+$", "\\1", x)) +
    scale_y_continuous(position = "right") +
    geom_text(data = distinct(eqtls_expression_df, locus, variant), 
              aes(x = 1.5, y = 3, label = variant), size = 3) +
    coord_cartesian(ylim = c(-3, 3.2)) +
    facet_wrap(~locus, ncol = 1, scales = "free") +
    labs(x = " ", y = " ") +
    theme_bw() +
    theme(text = element_text(size = 10, family = "Arial"),
          strip.text = element_blank())


grid1 <- plot_grid(legend, NULL, p1, p2, ncol = 2, 
                   rel_widths = c(2.5, 1), rel_heights = c(.07, 1))


tiff("./plots/Fig6.tiff", width = 6, height = 8, units = "in", res = 300)
ggdraw(grid1) + 
    draw_label("HLA lineage", 0.45, 0.025, size = 10, fontfamily = "Arial") +
    draw_label("eQTL genotype", 0.82, 0.025, size = 10, fontfamily = "Arial") +
    draw_label("PCA-corrected expression", .985, 0.5, size = 10, angle = 90, fontfamily = "Arial")
dev.off()

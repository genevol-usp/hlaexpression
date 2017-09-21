devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(ggpmisc)
library(ggplot2)

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
    select(gene_id, gene_name)

# genotype PCA 
make_pca_plot <- function(PC_x, PC_y) {
    
    ggplot(pcs, aes_string(PC_x, PC_y)) +
        geom_point(aes(color = pop)) +
        ggthemes::scale_color_colorblind() +
        theme_bw()
}

geuvadis_pops <- select(geuvadis_info, subject = name, pop)

pcs <- 
    read_delim("./pca_genotypes/eur.pca", delim = " ") %>%
    mutate(SampleID = sub("^.+_(PC\\d+)$", "\\1", SampleID)) %>%
    filter(SampleID %in% paste0("PC", 1:100)) %>%
    gather(subject, value, -1) %>%
    spread(SampleID, value) %>%
    mutate_at(vars(-subject), function(x) x/sqrt(sum(x^2))) %>%
    inner_join(geuvadis_pops, by = "subject")

png("./plots/genotype_pca.png", width = 10, height = 5, units = "in", res = 300)
p1 <- make_pca_plot("PC1", "PC2")
p2 <- make_pca_plot("PC2", "PC3")
p3 <- make_pca_plot("PC3", "PC4")
p4 <- make_pca_plot("PC4", "PC5")

leg <- get_legend(p1)

plot_grid(p1 + guides(color=FALSE), 
          p2 + guides(color=FALSE), 
          leg,
          p3 + guides(color=FALSE), 
          p4 + guides(color=FALSE), 
          nrow = 2, ncol = 3, rel_widths = c(3, 3, 1))
dev.off()

# Number of eQTLs according to index
pcs <- c(seq(0, 20, 5), seq(30, 100, 10))

egenes_pca_star_imgt <- 
    sprintf("./qtls_star/imgt/2-permutations/results/permutations_%d.significant.txt", pcs) %>%
    setNames(pcs) %>%
    parallel::mclapply(function(x) read_delim(x, delim = " ", col_names = FALSE), 
                       mc.cores = length(pcs)) %>%
    bind_rows(.id = "f") %>%
    count(f)

egenes_pca_star_pri <- 
    sprintf("./qtls_star/pri/2-permutations/results/permutations_%d.significant.txt", pcs) %>%
    setNames(pcs) %>%
    parallel::mclapply(function(x) read_delim(x, delim = " ", col_names = FALSE), 
                       mc.cores = length(pcs)) %>%
    bind_rows(.id = "f") %>%
    count(f)

egenes_df <- 
    list(imgt = egenes_pca_star_imgt, pri = egenes_pca_star_pri) %>%
    bind_rows(.id = "index") %>%
    mutate(f = as.integer(f)) %>%
    arrange(f, index)

png("./plots/n_of_egenes.png", width = 8, height = 5, units = "in", res = 300)
ggplot(egenes_df, aes(f, n, color = index, group = index)) + 
    geom_point(size = 2.5) + 
    geom_line() +
    scale_x_continuous(breaks = sort(unique(egenes_df$f))) +
    ggsci::scale_color_aaas() +
    theme_bw() +
    labs(x = "Number of PCs/factors", y = "Number of eGenes")
dev.off()

# Lineage and effects plot
## lineages
concordant <- 
    bind_rows(
        read_tsv("../expression/star/phase_hla_alleles/data/concordant_haps_classI.tsv") %>%
            gather(locus, allele, A:C), 
        read_tsv("../expression/star/phase_hla_alleles/data/concordant_haps_classII.tsv") %>%
            gather(locus, allele, DQA1:DRB1)) %>%
    arrange(subject, locus, hap)
    
haps_expression <-
    "../expression/star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv" %>%
    read_tsv() %>%
    filter(rank == 0)

hap_hla_genot <-
    haps_expression %>%
    select(subject, locus, hap, hla_allele)

hap_snps <-
    haps_expression %>%
    select(subject, locus, hap, variant_allele)

phen_best <-
    "./qtls_star/imgt/1-phenotypes/phenotypes_eur_60.bed.gz" %>%
    read_tsv() %>%
    inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
    select(gene_name, HG00096:NA20828) %>%
    mutate(gene_name = sub("HLA-", "", gene_name)) %>%
    gather(subject, resid, -gene_name) %>%
    select(subject, locus = gene_name, resid)

qtls_high_low <-
    left_join(hap_snps, phen_best, by = c("subject", "locus")) %>%
    inner_join(select(concordant, -allele)) %>%
    group_by(locus, variant_allele) %>% 
    summarize(eQTL = mean(resid)) %>%
    group_by(locus) %>%
    mutate(eQTL = ifelse(eQTL == max(eQTL), "High", "Low")) %>%
    ungroup()

lineage_df <-
    haps_expression %>%
    select(subject, locus, hla_allele, hap, tpm, variant_allele) %>%
    left_join(qtls_high_low, by = c("locus", "variant_allele")) %>%
    left_join(concordant, by = c("subject", "locus", "hap")) %>%
    mutate(eQTL = ifelse(is.na(allele), NA_character_, eQTL)) %>%
    select(subject, locus, hla_allele, eQTL, tpm) %>%
    mutate(hla_allele = hla_trimnames(hla_allele, 1),
           eQTL = factor(eQTL, levels = c("High", "Low")))

## effects
best_eqtl_locus <-
    read_tsv("../expression/star/phase_hla_alleles/best_eqtls.tsv") %>%
    filter(rank == 0L) %>%
    left_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(locus = gene_name, variant = var_id)

eqtl_info <- 
    read_tsv("../expression/star/phase_hla_alleles/best_eqtl_snps.vcf", comment = "##") %>%
    select(-`#CHROM`, -POS, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    gather(subject, genotype, -(ID:ALT)) %>%
    inner_join(best_eqtl_locus, by = c("ID" = "variant")) %>%
    separate(genotype, c("h1", "h2"), convert = TRUE) %>%
    gather(hap, allele, h1:h2) %>%
    mutate(locus = sub("HLA-", "", locus),
           allele = ifelse(allele == 0, REF, ALT)) %>%
    arrange(subject, locus, ID, allele) %>%
    group_by(subject, locus, ID) %>%
    summarize(genotype = paste(sort(allele), collapse = "/")) %>%
    select(subject, locus, variant = ID, genotype) %>%
    ungroup() %>%
    unite(id, genotype, locus, sep = "_", remove = FALSE)

eqtls_expression_df <- left_join(eqtl_info, phen_best, by = c("subject", "locus"))

## plot
png("./plots/lineage_and_effects.png", width = 12, height = 12, units = "in", res = 300)
p1 <- 
    ggplot(lineage_df, 
	   aes(x = reorder(hla_allele, tpm, FUN = median, na.rm = TRUE), 
	       y = tpm)) +
    geom_jitter(aes(color = eQTL), alpha = .5) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "grey5", alpha = .1) +
    scale_color_manual(values = c("Low" = "royalblue", "High" = "red"), 
		       na.value = "grey") +
    scale_x_discrete(labels = function(x) gsub("\\*", "*\n", x)) +
    facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
    labs(x = " ", y = "TPM") +
    theme_bw() +
    theme(axis.title = element_text(size = rel(1.2)),
	  axis.text = element_text(size = rel(1)),
	  legend.text = element_text(size = rel(1)),
	  strip.text = element_text(face = "bold", size = rel(1.2)),
	  legend.position = "top") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

p2 <-
    ggplot(eqtls_expression_df, aes(reorder(id, resid, "mean"), resid)) +
    geom_jitter(width = .25, alpha = 1/2, size = .75) +
    geom_smooth(aes(group = 1), method = lm, se = FALSE) +
    scale_x_discrete(labels = function(x) paste0(sub("^([^_]+).+$", "\\1", x), "\n")) +
    scale_y_continuous(position = "right") +
    geom_text(data = distinct(eqtls_expression_df, locus, variant), 
	      aes(x = 1.5, y = 3, label = variant)) +
    coord_cartesian(ylim = c(-3, 3.2)) +
    facet_wrap(~locus, ncol = 1, scales = "free") +
    labs(x = " ", y = " ") +
    theme_bw() +
    theme(axis.text = element_text(size = rel(1)),
	  strip.text = element_blank())

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

grid1 <- plot_grid(legend, NULL, p1, p2, ncol = 2, 
                   rel_widths = c(4, 1), rel_heights = c(.07, 1))

ggdraw(grid1) + 
    draw_label("HLA lineage", 0.44, 0.02, size = 16) +
    draw_label("1000G genotype", 0.88, 0.02, size = 16) +
    draw_label("PCA-corrected expession", .985, 0.5, size = 16, angle = 90)

dev.off()

# eQTL landscape around TSS
read_conditional <- function(path) {
      read_qtltools(path) %>%
      inner_join(select(gencode_hla, gene_id, gene_name), by = c("phen_id" = "gene_id")) %>%
      mutate(dist_tss = ifelse(strand == "+", 
			       var_from - phen_from,
			       phen_to - var_from),
	     nom_pval = -log10(bwd_pval)) %>%
      select(phen_id = gene_name, rank, var_id, var_from, var_to, dist_tss, nom_pval, 
	     slope = bwd_slope, best = bwd_best, signif = bwd_signif) %>%
      group_by(phen_id, var_id, best) %>%
      filter(rank == min(rank)) %>%
      ungroup() %>%
      mutate(rank = factor(rank))
}

plot_qtls <- function(conditional_df) {

  ggplot(conditional_df) +
    geom_vline(xintercept = 0, color = "grey25", size = 2) + 
    geom_point(data = filter(conditional_df, signif == 0), 
	       aes(dist_tss, nom_pval),
	       color = "grey", alpha = .1, show.legend = FALSE) +
    geom_point(data = filter(conditional_df, signif == 1L), 
	       aes(dist_tss, nom_pval, color = rank), 
	       alpha = .5) +
    geom_point(data = filter(conditional_df, best == 1L), 
	       aes(dist_tss, nom_pval, color = rank), 
	       size = 2) +
    geom_point(data = filter(conditional_df, best == 1L), 
	       aes(dist_tss, nom_pval), 
	       shape = 1, size = 2, color = "black", stroke = 1.5) +
    geom_hline(data = conditional_df %>% 
		 group_by(phen_id) %>% 
		 filter(signif == 0) %>% 
		 summarise(thres = max(nom_pval)),
	       aes(yintercept = thres), color = "black") +
    coord_cartesian(xlim = c(-1e6, +1e6)) +
    ggsci::scale_color_aaas() +
    scale_x_continuous(labels = scales::comma) +
    theme_minimal() +
    facet_wrap(~phen_id, scales = "free_y", ncol = 1) +
    labs(x = "distance from TSS", 
	 y = expression(paste("-log"[10], italic(Pvalue)))) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
}

conditional_star_imgt <-
    "./qtls_star/imgt/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional()

conditional_star_imgt %>%
    filter(best == 1L) %>%
    select(phen_id, rank, var_id, var_from, var_to, dist_tss, slope, nom_pval) %>%
    write_tsv("./plots/eqtls.tsv")

conditional_star_pri <-
    "./qtls_star/pri/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional()

png("./plots/qtls_landscape_imgt.png", height = 12, width = 10, units = "in", res = 300)
plot_qtls(conditional_star_imgt)
dev.off()

png("./plots/qtls_landscape_pri.png", height = 12, width = 10, units = "in", res = 300)
plot_qtls(conditional_star_pri)
dev.off()


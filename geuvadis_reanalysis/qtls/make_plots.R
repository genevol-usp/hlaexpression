devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(ggpmisc)
library(ggrepel)

hla_genes <- gencode_hla$gene_name

# genotype PCA 
make_pca_plot <- function(PC_x, PC_y) {
    
    ggplot(pcs, aes_string(PC_x, PC_y)) +
    geom_point(aes(color = pop)) +
	scale_color_manual(values = c("CEU" = "#8491B4", "FIN" = "#DC0000B2",
	                              "GBR" = "#7E6148FF", "TSI" = "#91D1C2")) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    theme(legend.text = element_text(size = 12)) +
    labs(color = "Population")
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
    sprintf("./star/supplemented/2-permutations/results/permutations_%d.significant.txt", pcs) %>%
    setNames(pcs) %>%
    plyr::ldply(. %>% read_delim(delim = " ", col_names = FALSE), .id = "f") %>%
    count(f)

egenes_pca_star_pri <- 
    sprintf("./star/transcriptome/2-permutations/results/permutations_%d.significant.txt", pcs) %>%
    setNames(pcs) %>%
    plyr::ldply(. %>% read_delim(, delim = " ", col_names = FALSE), .id = "f") %>%
    count(f)

egenes_df <- 
    list(imgt = egenes_pca_star_imgt, pri = egenes_pca_star_pri) %>%
    bind_rows(.id = "index") %>%
    mutate(f = as.integer(as.character(f))) %>%
    arrange(f, index)

png("./plots/n_of_egenes.png", width = 6, height = 3, units = "in", res = 300)
ggplot(egenes_df, aes(f, n, color = index, group = index)) + 
    geom_point(size = 2.5) + 
    geom_line() +
    scale_x_continuous(breaks = sort(unique(egenes_df$f))) +
    scale_color_manual(values = c(imgt = "#8491B4B2", pri = "#DC0000B2"),
		       labels = c(imgt = "HLA-personalized", pri = "Ref transcriptome")) +
    theme_bw() +
    labs(x = "Number of PCs", y = "Number of eGenes")
dev.off()

# eQTL landscape around TSS
read_conditional_hla <- function(path) {
    read_qtltools(path) %>%
        inner_join(hla_info_df, by = c("phen_id" = "gene_id")) %>%
        mutate(tss = ifelse(strand == "+", phen_from, phen_to),
               dist_to_tss = var_from - tss,
               pval = -log10(bwd_pval)) %>%
        select(gene = gene_name, strand, rank, var_id, var_from, dist, 
               dist_to_tss, pval, best = bwd_best, signif = bwd_signif) %>%
        group_by(gene, var_id) %>%
        filter(pval == max(pval)) %>%
        ungroup() %>%
        mutate(rank = factor(rank))
}

plot_qtls_index <- function(qtl_df) {
    
    tss_df <- qtl_df %>%
        distinct(index, gene, strand) %>%
        mutate(xend = ifelse(strand == "+", 5e4L, -5e4L)) %>%
        select(index, gene, xend)
    
    best_df <- filter(qtl_df, best == 1L)
    
    ggplot(data = qtl_df, aes(dist_to_tss, pval)) +
        geom_blank(data = qtl_df %>% 
                       group_by(gene, index) %>%
                       slice(which.max(pval)) %>% 
                       mutate(max_pval = pval + 1L),
                   aes(y = max_pval)) +
        coord_cartesian(xlim = c(-7e5, +7e5)) +
        geom_vline(aes(xintercept = 0), color = "grey", size = 1.5) + 
        geom_point(data = filter(qtl_df, signif == 0), 
                   aes(dist_to_tss, pval),
                   color = "grey", alpha = .1, show.legend = FALSE) +
        geom_segment(data = tss_df,
                     aes(x = 0, xend = xend, y = -1.5, yend = -1.5),
                     arrow = arrow(length = unit(0.2, "cm"), type = "closed", ends = "last"),
                     size = 1) +
        geom_point(data = filter(qtl_df, signif == 1L), 
                   aes(dist_to_tss, pval, color = rank), 
                   alpha = .25) +
        geom_point(data = best_df, 
                   aes(dist_to_tss, pval, color = rank), size = 2.5) +
        geom_point(data = best_df, 
                   aes(dist_to_tss, pval), 
                   shape = 1, size = 2.5, color = "black", stroke = 1) +
        geom_text_repel(data = best_df,
                        aes(dist_to_tss, pval, label = var_id),
                        fontface = "bold", size = 3.5, show.legend = FALSE) +
        scale_color_manual(values = c("0" = "#8491B4B2",
                                      "1" = "#DC0000B2",
                                      "2" = "gold3",
                                      "3" = "#7E6148FF",
                                      "4" = "#009E73")) +
        scale_x_continuous(breaks = seq(-5e5L, 5e5L, by = 2.5e5L),
                           labels = function(x) x/1e6L) +
        theme_bw() +
        theme(strip.text = element_text(size = 10, face = "bold"),
              strip.background = element_rect(fill = "grey"),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10)) +
        facet_grid(gene~index, scales = "free_y") +
        labs(x = "distance from TSS (Mb)", 
             y = expression(paste("-log"[10], italic(Pvalue)))) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
}

hla_info_df <- gencode_hla %>%
    select(gene_id, gene_name) 

hla_qtls <-
    "./star/supplemented/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_hla()

transcriptome_qtls <-
    "./star/transcriptome/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_hla()

hla_qtls %>%
    filter(best == 1L) %>%
    select(gene, rank, var_id, var_from, dist, pval) %>%
    write_tsv("./plots/eqtl.tsv")

png("./plots/qtls_landscape.png", height = 10, width = 10, units = "in", res = 300)
list("HLA-personalized" = hla_qtls, "Ref Transcriptome" = transcriptome_qtls) %>%
    bind_rows(.id = "index") %>%
    plot_qtls_index()
dev.off()


all_rank0 <- 
    "./star/supplemented/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(bwd_signif == 1L) %>%
    mutate(genome_context = ifelse(phen_id %in% gencode_hla$gene_id, "HLA", "genomewide")) %>%
    select(genome_context, dist, bwd_best)

png("./plots/qtls_density_geneStart.png", height = 2, width = 5, units = "in", res = 300)
ggplot(all_rank0) +
    coord_cartesian(xlim = c(-1e6, 1e6)) +
    geom_density(aes(x = dist, fill = genome_context), alpha = 1/2) +
    scale_fill_manual(values = c("genomewide" = "#8491B4B2", "HLA" = "#DC0000B2")) +
    theme_bw() +
    labs(x = "Gene start", fill = "Genome context")
dev.off()

# proportion of lead eQTLs further than the HLA-B lead eQTL
all_rank0 %>% 
    filter(bwd_best == 1L) %>%
    summarize(mean(abs(dist) >= 214757))


# Lineage and effects plot
## lineages
dist_to_ref <- "../../imgt_index_v2/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)
    
haps_expression <-
    "../expression/star/phase_hla_alleles/data/1000G_haps_expression_rsid.tsv" %>%
    read_tsv()

hap_hla_genot <- haps_expression %>%
    select(subject, locus, hap, hla_allele)

hap_snps <- haps_expression %>%
    select(subject, locus, hap, allele)

#calls <- haps_expression %>%
#    select(subject, locus, allele = hla_allele) %>%
#    mutate(locus = sub("HLA-", "", locus),
#           allele = hla_trimnames(allele, 1))
#
#gold <- pag %>%
#    mutate(allele = hla_trimnames(allele, 1))
#
#typing_errors <- calc_genotyping_accuracy(calls, gold, by_locus = FALSE) %>%
#    group_by(subject, locus) %>%
#    filter(any(!correct)) %>%
#    ungroup() %>%
#   distinct(subject, locus)

phen_best <- 
    "./star/supplemented/1-phenotypes/phenotypes_eur_60.bed.gz" %>%
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

lineage_df <- haps_expression %>%
    left_join(qtls_high_low, by = c("locus", "allele")) %>%
    mutate(eQTL = ifelse(grepl("/", hla_allele), NA, eQTL)) %>%
    select(subject, locus, hla_allele, eQTL, tpm) 

lineage_df_fixed <- lineage_df %>%
    filter(grepl("/", hla_allele)) %>%
    separate_rows(hla_allele, tpm, sep = "/") %>%
    distinct(subject, locus, hla_allele, eQTL, tpm) %>%
    bind_rows(filter(lineage_df, !grepl("/", hla_allele))) %>%
    arrange(subject, locus, hla_allele) %>%
    mutate(tpm = as.numeric(tpm)) %>%
    left_join(dist_to_ref, by = c("hla_allele" = "allele")) %>%
    mutate(lineage = hla_trimnames(hla_allele, 1),
           eQTL = factor(eQTL, levels = c("Low", "High")))

    
lineage_df_fixed %>%
    select(subject, locus, lineage, tpm) %>%
    group_by(lineage) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    split(.$locus) %>%
    map(~oneway.test(tpm~lineage, data = .) %>% broom::tidy()) %>%
    bind_rows(.id = "locus") %>%
    select(locus, num.df, denom.df, F = statistic, p.value) %>%
    write_tsv("./f_onewaytest_lineages.tsv")

lineage_df_fixed %>%
    select(subject, locus, lineage, tpm) %>%
    split(.$locus) %>%
    map(~lm(tpm~lineage, data = .) %>% anova() %>% broom::tidy()) %>%
    bind_rows(.id = "locus") %>%
    filter(term == "lineage") %>%
    select(locus, df, F = statistic, p.value) %>%
    write_tsv("./f_test_lineages.tsv")

## effects
best_eqtl_locus <- 
    read_tsv("./plots/eqtl.tsv") %>%
    filter(rank == 0) %>%
    select(locus = gene, variant = var_id)

eqtl_info <- 
    "../expression/star/phase_hla_alleles/best_eqtl_snps.vcf" %>%
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

## plot
png("./plots/lineage_and_effects.png", width = 12, height = 12, units = "in", res = 300)
p1 <- 
    ggplot(lineage_df_fixed, 
           aes(x = reorder(lineage, tpm, FUN = median, na.rm = TRUE), 
               y = tpm)) +
    geom_jitter(aes(color = eQTL)) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "grey5", alpha = .1) +
    scale_color_manual(values = c("Low" = "#8491B4", "High" = "#DC0000B2"), 
                       na.value = "grey") +
    scale_x_discrete(labels = function(x) sub("^(.+\\*)", "", x)) +
    facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
    labs(x = " ", y = "TPM", 
         color = "expression level associated with eQTL allele on same haplotype") +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          strip.text = element_text(face = "bold", size = 12),
          legend.position = "top") +
    guides(color = guide_legend(override.aes = list(size = 6)))

p2 <-
    ggplot(eqtls_expression_df, aes(reorder(id, resid, "mean"), resid)) +
    geom_jitter(width = .25, alpha = 1/2, size = .75) +
    geom_smooth(aes(group = 1), method = lm, se = FALSE) +
    scale_x_discrete(labels = function(x) sub("^([^_]+).+$", "\\1", x)) +
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
    draw_label("HLA lineage", 0.44, 0.01, size = 16) +
    draw_label("1000G genotype", 0.88, 0.01, size = 16) +
    draw_label("PCA-corrected expression", .985, 0.5, size = 16, angle = 90)
dev.off()


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
    ggsci::scale_color_npg() +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    theme(legend.text = element_text(size = 12)) +
    labs(color = "Population")
}

geuvadis_pops <- select(geuvadis_info, subject = name, pop)

pcs <- 
    read_delim("../data/pca_genotypes/eur.pca", delim = " ") %>%
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

pipelines <- c("Ref Transcriptome", "HLA-personalized",
               "Ref Transcriptome (quasi)", "HLA-personalized (quasi)", 
               "Conventional", "Ref Genome (Unique)") 

cols <- ggsci::pal_npg()(6) %>%
    setNames(pipelines)

pcs <- c(seq(0, 20, 5), seq(30, 100, 10))

egenes_hla_mapping <- 
    sprintf("./transcriptomemapping/hla_personalized/2-permutations/results/permutations_%d.significant.txt", pcs) %>%
    setNames(pcs) %>%
    map_df(~read_delim(., delim = " ", col_names = FALSE), .id = "f") %>%
    count(f)

egenes_ref_mapping <- 
    sprintf("./transcriptomemapping/reference/2-permutations/results/permutations_%d.significant.txt", pcs) %>%
    setNames(pcs) %>%
    map_df(~read_delim(., delim = " ", col_names = FALSE), .id = "f") %>%
    count(f)

egenes_df <- 
    list("HLA-personalized" = egenes_hla_mapping,
         "Ref Transcriptome" = egenes_ref_mapping) %>%
    bind_rows(.id = "index") %>%
    mutate(f = as.integer(as.character(f))) %>%
    arrange(f, index)

png("./plots/n_of_egenes.png", width = 6, height = 3, units = "in", res = 300)
ggplot(egenes_df, aes(f, n, color = index, group = index)) + 
    geom_point(size = 2.5) + 
    geom_line() +
    scale_x_continuous(breaks = sort(unique(egenes_df$f))) +
    scale_color_manual(values = cols) +
    theme_bw() +
    labs(x = "Number of PCs", y = "Number of eGenes")
dev.off()

# eQTL landscape around TSS
read_conditional_mhc <- function(file, mhc_genes) {
    qtltools_res <- 
        read_qtltools(file) %>%
        inner_join(mhc_genes, by = c("phen_id" = "gene_id")) %>%
        mutate(tss = ifelse(strand == "+", phen_from, phen_to),
               dist_to_tss = var_from - tss,
               pval = -log10(bwd_pval)) %>%
        select(gene = gene_name, tss, strand, rank, var_id, 
               var_from, dist, dist_to_tss, pval, best = bwd_best, 
               signif = bwd_signif) 
    
    fix_rank <- qtltools_res %>%
        filter(best == 1) %>%
        arrange(gene, desc(pval), rank) %>%
        group_by(gene) %>%
        mutate(new_rank = 1:n() - 1L) %>%
        ungroup() %>%
        select(gene, rank, new_rank)
    
    out_df <- left_join(qtltools_res, fix_rank, by = c("gene", "rank")) %>%
        select(gene, tss, strand, rank = new_rank, var_id, var_from, dist,
               dist_to_tss, pval, best, signif) %>%
        group_by(gene, var_id) %>%
        filter(pval == max(pval)) %>%
        ungroup() %>%
        mutate(gene = reorder(gene, tss),
            rank = as.character(rank))
    
    out_df
}

plot_qtls_index <- function(qtl_df) {
    
    tss_df <- qtl_df %>%
        distinct(index, gene, strand) %>%
        mutate(xend = ifelse(strand == "+", 1e5L, -1e5L)) %>%
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
                     aes(x = 0, xend = xend, y = -2.5, yend = -2.5),
                     arrow = arrow(length = unit(0.25, "cm"), type = "closed", ends = "last"),
                     size = 1, color = "darkgreen", alpha = .5) +
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
                        fontface = "bold", size = 4, segment.color = "gray35",
                        nudge_x = 1.5e5, force = 1,
                        show.legend = FALSE) +
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
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12)) +
        facet_grid(gene~index, scales = "free_y") +
        labs(x = "distance from TSS (Mb)", 
             y = expression(paste("-log"[10], italic(Pvalue)))) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
}

mhc_coords <- gencode_hla %>%
    summarise(start = min(start) - 1e6L, end = max(end) + 1e6L)

mhc_genes <- gencode_pri_gene %>%
    filter(start >= mhc_coords$start, end <= mhc_coords$end) %>%
    select(gene_id, gene_name)

hla_qtls <-
    "./transcriptomemapping/hla_personalized/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_mhc(mhc_genes)

ref_qtls <-
    "./transcriptomemapping/reference/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_mhc(mhc_genes)

qtls_df <- 
    list("HLA-personalized" = hla_qtls, "Ref Transcriptome" = ref_qtls) %>%
    bind_rows(.id = "index") %>%
    mutate(gene = reorder(gene, tss))

hla_qtls %>%
    filter(best == 1L, gene %in% hla_genes) %>%
    select(gene, rank, var_id, var_from, dist, pval) %>%
    write_tsv("./plots/eqtl.tsv")

hla_qtls <- hla_qtls %>%
    mutate(gene = reorder(gene, tss),
           colorize = ifelse(gene %in% gencode_hla$gene_name, 1L, 0L))

png("./plots/qtls_landscape.png", height = 12, width = 10, units = "in", res = 300)
plot_qtls_1 <- ggplot() +
    geom_point(data = filter(hla_qtls, colorize == 0L), 
               aes(var_from, pval), 
               color = "grey", alpha = .1, size = .5) +
    geom_point(data = filter(hla_qtls, colorize == 1L),
               aes(var_from, pval, color = gene),
               alpha = .25, size = .5) +
    ggsci::scale_color_npg() +
    scale_x_continuous(labels = function(x) x/1e6L) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    coord_cartesian(xlim = c(29.2e6L, 33.8e6L)) +
    labs(x = "position (Mb)",
         y = expression(paste("-log"[10], italic(Pvalue))))

plot_qtls_2 <- qtls_df %>%
    filter(gene %in% gencode_hla$gene_name) %>%
    plot_qtls_index()

plot_grid(plot_qtls_1, plot_qtls_2, ncol = 1,
          rel_heights = c(.25, 1),
          labels = c("A", "B"))
dev.off()


# QTLs density genome-wide vs HLA
all_rank0 <- 
    "./transcriptomemapping/hla_personalized/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(bwd_signif == 1L) %>%
    mutate(genome_context = ifelse(phen_id %in% gencode_hla$gene_id, "HLA", "genomewide")) %>%
    select(genome_context, dist, bwd_best)
#
#png("./plots/qtls_density_geneStart.png", height = 2, width = 5, units = "in", res = 300)
#ggplot(all_rank0) +
#    coord_cartesian(xlim = c(-1e6, 1e6)) +
#    geom_density(aes(x = dist, fill = genome_context), alpha = 1/2) +
#    scale_fill_manual(values = c("genomewide" = "#8491B4B2", "HLA" = "#DC0000B2")) +
#    theme_bw() +
#    labs(x = "Gene start", fill = "Genome context")
#dev.off()

# proportion of lead eQTLs further than the HLA-B lead eQTL
all_rank0 %>% 
    filter(bwd_best == 1L) %>%
    summarize(mean(abs(dist) >= 214757))


# Lineage and effects plot
## lineages
lineage_df <- 
    "../expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele, tpm) %>%
    mutate(subject = convert_ena_ids(subject),
           allele = gsub("IMGT_", "", allele),
           allele_3f = hla_trimnames(allele, 3),
           allele_2f = hla_trimnames(allele, 2),
           lineage = hla_trimnames(allele, 1))

lineage_df %>%
    select(subject, locus, lineage, tpm) %>%
    group_by(lineage) %>%
    filter(n() > 10) %>%
    ungroup() %>%
    split(.$locus) %>%
    map(~oneway.test(tpm~lineage, data = .) %>% broom::tidy()) %>%
    bind_rows(.id = "locus") %>%
    select(locus, num.df, denom.df, F = statistic, p.value) %>%
    write_tsv("./f_onewaytest_lineages.tsv")

lineage_df %>%
    select(subject, locus, lineage, tpm) %>%
    group_by(lineage) %>%
    filter(n() > 10) %>%
    ungroup() %>%
    split(.$locus) %>%
    map(~lm(tpm~lineage, data = .) %>% anova() %>% broom::tidy()) %>%
    bind_rows(.id = "locus") %>%
    filter(term == "lineage") %>%
    select(locus, df, F = statistic, p.value) %>%
    write_tsv("./f_test_lineages.tsv")

dist_to_ref <- "../../imgt_index_v2/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)

haps_expression <- "../phase_hla/phase_hla_haps_snps.tsv" %>%
    read_tsv() %>%
    filter(rank == 0L) %>%
    rename(hla_allele = allele_gene) %>%
    left_join(lineage_df, by = c("subject", "locus", "hla_allele" = "allele_3f")) %>%
    distinct() %>%
    group_by(subject, locus) %>%
    filter(n() != 4L) %>%
    ungroup()

hap_hla_genot <- haps_expression %>%
    select(subject, locus, hap, hla_allele)

hap_snps <- haps_expression %>%
    select(subject, locus, hap, allele = allele_snp)

phen_best <- 
    "./transcriptomemapping/hla_personalized/1-phenotypes/phenotypes_eur_60.bed.gz" %>%
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
    select(subject, locus, allele, eQTL, tpm) %>%
    left_join(dist_to_ref, by = c("allele")) %>%
    mutate(lineage = hla_trimnames(allele, 1),
           allele_2f = hla_trimnames(allele, 2),
           eQTL = factor(eQTL, levels = c("Low", "High"))) %>%
    group_by(lineage) %>%
    filter(n() >= 10L) %>%
    ungroup() %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name))
    

## effects
best_eqtl_locus <- read_tsv("./plots/eqtl.tsv") %>%
    filter(rank == 0) %>%
    select(locus = gene, variant = var_id)

eqtl_info <- "../phase_hla/eqtl_snps.vcf" %>%
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

## plot
plot_lineages <- function(df) { 
    ggplot(data = df,
           aes(x = reorder(lineage, tpm, FUN = median, na.rm = TRUE), 
               y = tpm)) +
    geom_jitter(aes(color = eQTL), size = .75) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "grey5", alpha = .1) +
    scale_color_manual(values = c("Low" = ggsci::pal_npg()(6)[6], 
                                  "High" = ggsci::pal_npg()(1), 
                                  "ND" = "grey")) +
    scale_x_discrete(labels = function(x) sub("^(.+\\*)", "", x)) +
    facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
    labs(x = " ", y = "TPM") +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          strip.text = element_text(face = "bold", size = 12),
          legend.position = "top") +
    guides(color = guide_legend(override.aes = list(size = 6)))
}
    
plot_slopes <- function(df) {
    ggplot(df, aes(reorder(id, resid, "mean"), resid)) +
    geom_jitter(width = .25, alpha = 1/2, size = .75) +
    geom_smooth(aes(group = 1), method = lm, se = FALSE) +
    scale_x_discrete(labels = function(x) sub("^([^_]+).+$", "\\1", x)) +
    scale_y_continuous(position = "right") +
    geom_text(data = distinct(df, locus, variant), 
              aes(x = 1.5, y = 3, label = variant)) +
    coord_cartesian(ylim = c(-3, 3.2)) +
    facet_wrap(~locus, ncol = 1, scales = "free") +
    labs(x = " ", y = " ") +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          strip.text = element_blank())
}

png("./plots/lineage_and_effects.png", width = 7, height = 10, units = "in", res = 300)
p1 <- plot_lineages(lineage_phased)

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- plot_slopes(eqtls_expression_df)

grid1 <- plot_grid(legend, NULL, p1, p2, ncol = 2, 
                   rel_widths = c(2.5, 1), rel_heights = c(.07, 1))

ggdraw(grid1) + 
    draw_label("HLA lineage", 0.45, 0.01, size = 14) +
    draw_label("eQTL genotype", 0.82, 0.01, size = 14) +
    draw_label("PCA-corrected expression", .985, 0.5, size = 14, angle = 90)
dev.off()

# CaVEMaN
caveman <- "./transcriptomemapping/hla_personalized/caveman/results.hla" %>%
    read_tsv() %>%
    mutate(gene = factor(gene, levels = gencode_hla$gene_name),
           rank = as.character(rank))

png("./plots/caveman.png", width = 10, height = 3, units = "in", res = 300)
ggplot(data = caveman, aes(rank, Probability, fill = index)) +
    geom_bar(stat = "identity", position = "dodge", alpha = .8) +
    ggsci::scale_fill_npg() +
    facet_wrap(~gene, nrow = 1, scales = "free_x") +
    labs(x = "") +
    theme(legend.position = "top")
dev.off()

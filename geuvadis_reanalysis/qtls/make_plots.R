devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(ggpmisc)

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
        inner_join(select(gencode_hla, gene_id, gene_name), 
		   by = c("phen_id" = "gene_id")) %>%
        mutate(dist_tss = ifelse(strand == "+", var_from - phen_from, phen_to - var_from),
               pval = -log10(bwd_pval)) %>%
        select(gene = gene_name, rank, var_id, var_from, dist, dist_tss, 
	       pval, slope = bwd_slope, best = bwd_best, signif = bwd_signif) %>%
        group_by(gene, var_id) %>%
        filter(pval == max(pval)) %>%
        ungroup() %>%
        mutate(rank = factor(rank))
}

plot_qtls <- function(conditional_df) {
    
    ggplot(conditional_df) +
        geom_vline(xintercept = 0, color = "grey45", size = 2) + 
        geom_point(data = filter(conditional_df, signif == 0), 
                   aes(dist_tss, pval),
                   color = "grey", alpha = .1, show.legend = FALSE) +
        geom_point(data = filter(conditional_df, signif == 1L), 
                   aes(dist_tss, pval, color = rank), 
                   alpha = .5) +
        geom_hline(data = conditional_df %>% 
                       group_by(gene) %>% 
                       filter(signif == 0) %>% 
                       summarise(thres = max(pval)),
                   aes(yintercept = thres), color = "black") +
        geom_point(data = filter(conditional_df, best == 1L), 
                   aes(dist_tss, pval, color = rank), size = 2.5) +
        geom_point(data = filter(conditional_df, best == 1L), 
                   aes(dist_tss, pval), 
                   shape = 1, size = 2.5, color = "black", stroke = 1) +
        ggrepel::geom_label_repel(data = filter(conditional_df, best == 1L) %>%
                                      mutate(lab = paste0(var_id, " (", Probability, ")")),
                                  aes(dist_tss, pval, label = lab, fill = rank),
                                  fontface = "bold", size = 3, show.legend = FALSE) +
        coord_cartesian(xlim = c(-1e6, +1e6)) +
        scale_color_manual(values = c("0" = "#8491B4B2",
                                      "1" = "#DC0000B2",
                                      "2" = "gold3",
                                      "3" = "#7E6148FF",
                                      "4" = "#009E73")) +
        scale_fill_manual(values = c("0" = "#8491B4B2",
                                     "1" = "#DC0000B2",
                                     "2" = "gold3",
                                     "3" = "#7E6148FF",
                                     "4" = "#009E73")) +
        scale_x_continuous(labels = function(x) x/1e6L) +
        theme_minimal() +
        theme(strip.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10)) +
        facet_wrap(~gene, scales = "free_y", ncol = 1) +
        geom_blank(data = conditional_df %>% 
                       group_by(gene) %>%
                       slice(which.max(pval)) %>% 
                       mutate(max_pval = pval + 2.5),
                   aes(y = max_pval)) +
        labs(x = "distance from TSS (Mb)", 
             y = expression(paste("-log"[10], italic(Pvalue)))) +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
}


plot_qtls_indices <- function(conditional_df) {

  ggplot(conditional_df) +
    geom_vline(xintercept = 0, color = "grey45", size = 2) + 
    geom_point(data = filter(conditional_df, signif == 0), 
               aes(dist_tss, pval),
               color = "grey", alpha = .1, show.legend = FALSE) +
    geom_point(data = filter(conditional_df, signif == 1L), 
               aes(dist_tss, pval, color = rank), 
               alpha = .5) +
    geom_hline(data = conditional_df %>% 
                   group_by(gene) %>% 
                   filter(signif == 0) %>% 
                   summarise(thres = max(pval)),
               aes(yintercept = thres), color = "black") +
    geom_point(data = filter(conditional_df, best == 1L), 
               aes(dist_tss, pval, color = rank), size = 2.5) +
    geom_point(data = filter(conditional_df, best == 1L), 
               aes(dist_tss, pval), 
               shape = 1, size = 2.5, color = "black", stroke = 1) +
    ggrepel::geom_label_repel(data = filter(conditional_df, best == 1L) %>%
                                  mutate(lab = paste0(var_id, " (", Probability, ")")),
                              aes(dist_tss, pval, label = lab, fill = rank),
                              fontface = "bold", size = 3, show.legend = FALSE) +
    coord_cartesian(xlim = c(-1e6, +1e6)) +
    scale_color_manual(values = c("0" = "#8491B4B2",
                                  "1" = "#DC0000B2",
                                  "2" = "#F0E442B2",
                                  "3" = "#7E6148B2",
                                  "4" = "#009E73")) +
    scale_fill_manual(values = c("0" = "#8491B4B2",
                                 "1" = "#DC0000B2",
                                 "2" = "#F0E442B2",
                                 "3" = "#7E6148B2",
                                 "4" = "#009E73")) +
    scale_x_continuous(labels = function(x) x/1e6L) +
    theme_minimal() +
    theme(strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)) +
    facet_grid(gene~index, scales = "free_y") +
    geom_blank(data = conditional_df %>% 
                   group_by(gene) %>%
                   slice(which.max(pval)) %>% 
                   mutate(max_pval = pval + 2.5),
               aes(y = max_pval)) +
    labs(x = "distance from TSS (Mb)", 
         y = expression(paste("-log"[10], italic(Pvalue)))) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
}

caveman_scores <- read_tsv("./star/supplemented/caveman/results.hla") %>%
    select(index, gene_name, rank, var_id, Probability) %>%
    mutate(rank = as.character(rank), Probability = round(Probability, digits = 2))

conditional_star_imgt <-
    "./star/supplemented/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_hla()

conditional_star_imgt %>%
    filter(best == 1L) %>%
    select(gene, rank, var_id, var_from, dist, dist_tss, slope, pval) %>%
    write_tsv("./plots/eqtl.tsv")

#conditional_star_imgt_v2 <-
#    "./star/supplemented_fix_genos/3-conditional_analysis/conditional_60_all.txt.gz" %>%
#    read_conditional_hla()

#"./star/supplemented/3-conditional_analysis/conditional_60_all.txt.gz" %>%
#    read_qtltools() %>%
#    filter(phen_id == gencode_hla$gene_id[gencode_hla$gene_name == "HLA-C"],
#	   rank == 0) %>%
#    filter(bwd_pval == min(bwd_pval))

#conditional_star_imgt_v2 %>%
#    filter(best == 1L) %>%
#    select(gene, rank, var_id, var_from, dist, dist_tss, slope, pval)

conditional_star_pri <-
    "./star/transcriptome/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_hla()

conditional_df <- 
    list("HLA_personalized" = conditional_star_imgt,
         "Reference" = conditional_star_pri) %>%
    bind_rows(.id = "index") %>%
    left_join(caveman_scores, 
              by = c("index", "gene" = "gene_name", "rank", "var_id")) %>%
    mutate(rank = factor(rank))

png("./plots/qtls_landscape_imgt.png", height = 10, width = 8, units = "in", res = 300)
conditional_df %>%
    filter(index == "HLA_personalized") %>%
    plot_qtls()
dev.off()

png("./plots/qtls_landscape_pri.png", height = 10, width = 10, units = "in", res = 300)
conditional_df %>%
    filter(index == "Reference") %>%
    plot_qtls()
dev.off()

png("./plots/qtls_landscape.png", height = 10, width = 12, units = "in", res = 300)
plot_qtls_indices(conditional_df)
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
dist_to_ref <- "../../simulation/PEreads_75bp/data/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)
    
haps_expression <-
    "../expression/star/phase_hla_alleles/data/1000G_haps_expression_rsid.tsv" %>%
    read_tsv()

hap_hla_genot <- haps_expression %>%
    select(subject, locus, hap, hla_allele)

hap_snps <- haps_expression %>%
    select(subject, locus, hap, allele)

phen_best <- "./star/supplemented/1-phenotypes/phenotypes_eur_60.bed.gz" %>%
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
best_eqtl_locus <- read_tsv("./plots/eqtl.tsv") %>%
    filter(rank == 0) %>%
    select(locus = gene, variant = var_id)

eqtl_info <- "../expression/star/phase_hla_alleles/best_eqtl_snps.vcf" %>%
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
    scale_color_manual(values = c("Low" = "#8491B4B2", "High" = "#DC0000B2"), 
                       na.value = "grey") +
    scale_x_discrete(labels = function(x) gsub("\\*", "*\n", x)) +
    facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
    labs(x = " ", y = "TPM") +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          strip.text = element_text(face = "bold", size = 12),
          legend.position = "top") +
    guides(color = guide_legend(override.aes = list(size = 4)))

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
    draw_label("PCA-corrected expression", .985, 0.5, size = 16, angle = 90)
dev.off()


devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(ggpmisc)

hla_genes <- sort(gencode_hla$gene_name)

# genotype PCA 
make_pca_plot <- function(PC_x, PC_y) {
    
    ggplot(pcs, aes_string(PC_x, PC_y)) +
        geom_point(aes(color = pop)) +
	scale_color_manual(values = c("CEU" = "#8491B4B2",
				      "FIN" = "#DC0000B2",
				      "GBR" = "#7E6148FF",
				      "TSI" = "#91D1C2B2")) +
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
    sprintf("./star/imgt/2-permutations/results/permutations_%d.significant.txt", pcs) %>%
    setNames(pcs) %>%
    parallel::mclapply(function(x) read_delim(x, delim = " ", col_names = FALSE), 
                       mc.cores = length(pcs)) %>%
    bind_rows(.id = "f") %>%
    count(f)

egenes_pca_star_pri <- 
    sprintf("./star/pri/2-permutations/results/permutations_%d.significant.txt", pcs) %>%
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

png("./plots/n_of_egenes.png", width = 6, height = 4, units = "in", res = 300)
ggplot(egenes_df, aes(f, n, color = index, group = index)) + 
    geom_point(size = 2.5) + 
    geom_line() +
    scale_x_continuous(breaks = sort(unique(egenes_df$f))) +
    ggsci::scale_color_aaas() +
    theme_bw() +
    labs(x = "Number of PCs", y = "Number of eGenes")
dev.off()

# eQTL landscape around TSS
read_conditional_hla <- function(path) {
    read_qtltools(path) %>%
        inner_join(select(gencode_hla, gene_id, gene_name, tss), 
		   by = c("phen_id" = "gene_id")) %>%
        mutate(dist_tss = ifelse(strand == "+", var_from - tss, tss - var_from),
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
        ggrepel::geom_label_repel(data = filter(conditional_df, best == 1L),
                                  aes(dist_tss, pval, label = var_id, fill = rank),
                                  fontface = 'bold', size = 3.5,  
                                  show.legend = FALSE) +
        coord_cartesian(xlim = c(-1e6, +1e6)) +
        scale_color_manual(values = c("0" = "#8491B4B2",
                                      "1" = "#DC0000B2",
                                      "2" = "gold3",
                                      "3" = "black",
                                      "4" = "#009E73")) +
        scale_fill_manual(values = c("0" = "#8491B4B2",
                                      "1" = "#DC0000B2",
                                      "2" = "gold3",
                                      "3" = "black",
                                      "4" = "#009E73")) +
        scale_x_continuous(labels = scales::comma) +
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
        labs(x = "distance from TSS", 
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
    ggrepel::geom_text_repel(data = filter(conditional_df, best == 1L),
                             aes(dist_tss, pval, label = Probability)) +
    coord_cartesian(xlim = c(-1e6, +1e6)) +
    scale_color_manual(values = c("0" = "#8491B4B2",
                                  "1" = "#DC0000B2",
                                  "2" = "#F0E442B2",
                                  "3" = "black",
                                  "4" = "#009E73")) +
    scale_x_continuous(labels = scales::comma) +
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
    labs(x = "distance from TSS", 
         y = expression(paste("-log"[10], italic(Pvalue)))) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
}

caveman_scores <- read_tsv("./star/imgt/caveman/results.hla") %>%
    select(index, gene_name, rank, var_id, Probability) %>%
    mutate(rank = as.character(rank), Probability = round(Probability, digits = 2))

conditional_star_imgt <-
    "./star/imgt/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_hla()

#"./star/imgt/3-conditional_analysis/conditional_60_all.txt.gz" %>%
#    read_qtltools() %>%
#    filter(phen_id == gencode_hla$gene_id[gencode_hla$gene_name == "HLA-C"],
#	   rank == 0) %>%
#    filter(bwd_pval == min(bwd_pval))

conditional_star_imgt %>%
    filter(best == 1L) %>%
    select(phen_id, rank, var_id, var_from, dist, dist_tss, slope, nom_pval) %>%
    write_tsv("./plots/eqtl.tsv")

conditional_star_pri <-
    "./star/pri/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_hla()

png("./plots/qtls_landscape_imgt.png", height = 10, width = 8, units = "in", res = 300)
plot_qtls(conditional_star_imgt)
dev.off()

png("./plots/qtls_landscape_pri.png", height = 10, width = 8, units = "in", res = 300)
plot_qtls(conditional_star_pri)
dev.off()

conditional_df <- 
    list("HLA_personalized" = conditional_star_imgt,
         "Reference" = conditional_star_pri) %>%
    bind_rows(.id = "index") %>%
    left_join(caveman_scores, 
              by = c("index", "gene" = "gene_name", "rank", "var_id")) %>%
    mutate(rank = factor(rank))

png("./plots/qtls_landscape.png", height = 10, width = 12, units = "in", res = 300)
plot_qtls_indices(conditional_df)
dev.off()


all_rank0 <- 
    "./star/imgt/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_qtltools() %>%
    filter(bwd_signif == 1L) %>%
    mutate(genome_context = ifelse(phen_id %in% gencode_hla$gene_id, "HLA", "genomewide")) %>%
    select(genome_context, dist, bwd_best)

png("./plots/qtls_density_geneStart.png", height = 3, width = 5, units = "in", res = 300)
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


## CRD
gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene_name, start, end, strand)

crd <- "../../LCL_ALL.chr6.subset.txt.gz" %>%
    read_delim(col_names = FALSE, delim = " ") %>%
    select(-X2, -X4, -X6, -X8, -X9, -X11) %>%
    filter(between(X3, min(gencode_hla_v19$start) - 1e6, max(gencode_hla_v19$end) + 1e6),
           between(X7, min(gencode_hla_v19$start) - 1e6, max(gencode_hla_v19$end) + 1e6)) %>%
    mutate(x = (X1 + X5)/2L,
           x2 = (X3 + X7)/2L,
           y = X5 - X1,
           r2 = X10^2) %>%
    select(-X10)

gene_pos <- gencode_hla_v19 %>%
    mutate(pos = ifelse(strand == "+", start, end),
           x1 = crd$x[map_dbl(pos, ~which.min(abs(. - crd$x2)))])

png("./plots/crd.png", width = 12, height = 5, units = "in", res = 200)
ggplot() +
    geom_point(data = crd, aes(x = x, y = y, alpha = r2), color = "blue") +
    scale_alpha_continuous(range = c(0, 1)) +
    scale_x_continuous(breaks = crd$x[seq(1, nrow(crd), 2e4)],
                       labels = round(crd$x2[seq(1, nrow(crd), 2e4)]/1e6L, digits = 1)) +
    geom_point(data = filter(gene_pos, strand == "+"), 
               aes(x = x1, y = -5), shape = ">", size = 5, color = "red") +
    geom_point(data = filter(gene_pos, strand == "-"), 
               aes(x = x1, y = -20), shape = "<", size = 5, color = "red") +
    geom_text(data = filter(gene_pos, strand == "+"), 
              aes(x = x1, y = -10, label = gene_name),
              size = 3) +
    ggrepel::geom_text_repel(data = filter(gene_pos, strand == "-"), 
                             aes(x = x1, y = -25, label = gene_name),
                             size = 3) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(color = "black")) +
    labs(x = "Genomic position (Mb)") +
    coord_cartesian(ylim = c(-25, max(crd$y)))
dev.off()


# Lineage and effects plot
## lineages
concordant <- 
    bind_rows(
        read_tsv("../expression/star/phase_hla_alleles/data/concordant_haps_classI.tsv") %>%
            gather(locus, allele, A:C), 
        read_tsv("../expression/star/phase_hla_alleles/data/concordant_haps_classII.tsv") %>%
            gather(locus, allele, DPB1:DRB1)) %>%
    arrange(subject, locus, hap)

dist_to_ref <- "../../simulation/PEreads_75bp/data/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)
    
haps_expression <-
    "../expression/star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv" %>%
    read_tsv()

hap_hla_genot <- haps_expression %>%
    select(subject, locus, hap, hla_allele)

hap_snps <- haps_expression %>%
    select(subject, locus, hap, variant_allele)

phen_best <- "./star/imgt/1-phenotypes/phenotypes_eur_60.bed.gz" %>%
    read_tsv() %>%
    inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
    select(gene_name, HG00096:NA20828) %>%
    mutate(gene_name = sub("HLA-", "", gene_name)) %>%
    gather(subject, resid, -gene_name) %>%
    select(subject, locus = gene_name, resid)

qtls_high_low <- left_join(hap_snps, phen_best, by = c("subject", "locus")) %>%
    inner_join(select(concordant, -allele)) %>%
    group_by(locus, variant_allele) %>% 
    summarize(eQTL = mean(resid)) %>%
    group_by(locus) %>%
    mutate(eQTL = ifelse(eQTL == max(eQTL), "High", "Low")) %>%
    ungroup()

lineage_df <- haps_expression %>%
    select(subject, locus, hla_allele, hap, tpm, variant_allele) %>%
    left_join(qtls_high_low, by = c("locus", "variant_allele")) %>%
    left_join(concordant, by = c("subject", "locus", "hap")) %>%
    mutate(eQTL = ifelse(is.na(allele), NA_character_, eQTL)) %>%
    select(subject, locus, hla_allele, eQTL, tpm) %>%
    left_join(dist_to_ref, by = c("hla_allele" = "allele")) %>%
    mutate(hla_allele = hla_trimnames(hla_allele, 1),
           eQTL = factor(eQTL, levels = c("Low", "High")))

lineage_exp <- haps_expression %>%
    mutate(lineage = hla_trimnames(hla_allele, 1)) %>%
    select(subject, locus, lineage, tpm) 

lineage_exp %>%
    group_by(lineage) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    split(.$locus) %>%
    map(~oneway.test(tpm~lineage, data = .) %>% broom::tidy()) %>%
    bind_rows(.id = "locus") %>%
    select(locus, num.df, denom.df, F = statistic, p.value) %>%
    write_tsv("./f_onewaytest_lineages.tsv")

lineage_exp %>%
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
    select(locus = phen_id, variant = var_id)

eqtl_info <- "../expression/star/phase_hla_alleles/best_eqtl_snps.vcf" %>%
    read_tsv(comment = "##") %>%
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


# mapping bias

#png("./plots/regulation_by_divergence.png", width = 10, height = 5, units = "in", res = 200)
#lineage_df %>%
#    filter(!is.na(eQTL)) %>%
#    ggplot() +
#    geom_boxplot(aes(eQTL, dist, color = eQTL), fill = NA) +
#    scale_color_manual(values = c("Low" = "#8491B4", "High" = "#DC0000")) +
#    facet_wrap(~locus) +
#    theme_bw() +
#    labs(x = "", y = "divergence to reference (%)")
#dev.off()
#

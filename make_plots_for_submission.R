devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(cowplot)
library(GGally)
library(scales)
library(ggpmisc)
library(ggrepel)

# Fig2

pipelines <- 
    c("HLApers (ML)", "Ref Transcriptome (ML)", "Ref Genome (Unique)") 

cols <- ggsci::pal_npg()(4) %>%
    .[c(1, 2, 4)] %>%
    setNames(pipelines)

allele_dist <- "./imgt_index/distance_to_ref/distances_to_reference.tsv" %>%
    read_tsv()

hla_genes <- gencode_hla$gene_name

hla_regex <- sub("HLA-", "", hla_genes) %>%
    paste(collapse = "|") %>%
    paste0("IMGT_(", ., ")")

ground_truth <- read_tsv("./simulation/data/phenotypes_trueCounts.tsv") %>%
    filter(grepl(hla_regex, Name)) %>%
    mutate(locus = sub("^IMGT_([^*]+).+$", "HLA-\\1", Name)) %>%
    group_by(subject, locus) %>%
    summarize(true_counts = sum(TrueCounts)) %>%
    ungroup()

uniq_reads <-
    "./simulation/expression/1-map_to_genome/quantifications_uniq/gw_quants.tsv" %>%
    read_tsv(col_names = FALSE) %>%
    inner_join(gencode_hla, c("X2" = "gene_id")) %>%
    select(subject = X1, locus = gene_name, est_counts = X3)

transcriptMap_ref <- 
    "./simulation/expression/3-map_to_transcriptome/reference/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv()

transcriptMap_pers <- 
    "./simulation/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    mutate(allele = gsub("IMGT_", "", allele)) %>%
    left_join(allele_dist, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), 
              tpm = sum(tpm), 
              dist = mean(dist)) %>%
    ungroup()

hla_counts_df <-
    left_join(transcriptMap_pers, transcriptMap_ref, by = c("subject", "locus")) %>%
    left_join(uniq_reads, by = c("subject", "locus")) %>%
    select(subject, locus, 
           counts.transcriptMap_pers = est_counts.x, 
           counts.transcriptMap_ref = est_counts.y,
           counts.uniq = est_counts,
           dist) %>%
    left_join(ground_truth, by = c("subject", "locus")) %>%
    gather(pipeline, counts, starts_with("counts")) %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name),
           pipeline = sub("^counts\\.", "", pipeline),
           pipeline = recode(pipeline, 
                             transcriptMap_pers = "HLApers (ML)",
                             transcriptMap_ref = "Ref Transcriptome (ML)",
                             uniq = "Ref Genome (Unique)"),
           pipeline = factor(pipeline, levels = pipelines),
           prop_mapped = counts/true_counts)

tiff("./plots_for_submission/Fig2.tiff", width = 5, height = 4, units = "in", res = 300)
ggplot(hla_counts_df, aes(dist, prop_mapped, color = pipeline)) +
    geom_point(size = 1) +
    geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
              alpha = 1/2, size = 1) +
    scale_x_continuous(labels = percent,
                       breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    facet_wrap(~locus, scales = "free_x") +
    scale_color_manual(values = cols) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Arial"),
          strip.text = element_text(face = "bold"),
          panel.spacing.x = unit(1, "lines"),
          legend.position = "top") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(x = "Sequence divergence to the HLA reference allele", 
         y = expression(frac(Aligned~reads, Simulated~reads)))
dev.off()


# Fig3
 
hla_genes <- gencode_hla$gene_name

geuvadis_ids <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(subject = ena_id, name)

imgt_loci <- readLines("./imgt_index/imgt_loci.txt") %>%
    paste0("HLA-", .)

mapping_imgt <- 
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% imgt_loci) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    group_by(locus) %>%
    filter(mean(tpm) >= 100) %>%
    ungroup() %>%
    mutate(locus = reorder(locus, tpm, median),
           locus = factor(locus, levels = rev(levels(locus))))

tiff("./plots_for_submission/Fig3.tiff", width = 5, height = 3, units = "in", res = 300)
ggplot(mapping_imgt, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = comma) +
    theme_bw() +
    theme(text = element_text(size = 9, family = "Arial"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = "TPM")
dev.off()


# Fig4
 
scatter_plot_cors <- function(df, x_var, y_var, dist_var) {
    
    cor_df <- df %>%
        group_by(locus) %>%
        do(data.frame(r = cor(.[[x_var]], .[[y_var]]),
                      p = cor(.[[x_var]], .[[y_var]], method = "spearman"),
                      x = min(.[[x_var]]),
                      y = max(.[[y_var]]))) %>%
        ungroup() %>%
        mutate_at(vars(r, p), ~round(., digits = 2)) %>%
        mutate(label = paste("r == ", r))

    ggplot(df, aes_string(x_var, y_var)) +
        geom_abline() +
        geom_point(aes_string(color = dist_var), size = .8) +
        scale_x_continuous(labels = comma, breaks = pretty_breaks(2)) +
        scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) +
        scale_color_gradient(low = "white", high = "darkred") +
        geom_text(data = cor_df, aes(x, y, label = label), parse = TRUE, 
                  hjust = "inward", vjust = "inward", size = 3, 
                  color = "white") +
        facet_wrap(~locus, scales = "free") +
        theme_dark() +
        theme(text = element_text(size = 9, family = "Arial"),
              strip.text = element_text(face = "bold"),
              legend.position = c(.75, .13)) +
        guides(color = guide_colorbar(barheight = 4))
}


mapping_ref <-
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/reference/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    select(subject, locus, tpm)

imgt_df <- list(mapping.hla_personalized = mapping_imgt, 
                mapping.ref = mapping_ref) %>%
    bind_rows(.id = "pipeline")

dist_ref <- read_tsv("./imgt_index/distances_to_reference.tsv")

mapping_imgt_dist <- 
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    left_join(dist_ref, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarise(dist = mean(dist)) %>%
    ungroup()

imgt_df_wide <- imgt_df %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    spread(pipeline, tpm) %>%
    left_join(mapping_imgt_dist) %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name))

tiff("./plots_for_submission/Fig4.tiff", width = 5.2, height = 3.5, units = "in", res = 300)
scatter_plot_cors(imgt_df_wide, "mapping.hla_personalized", "mapping.ref", "dist") +
    labs(x = "HLApers", y = "Ref Transcriptome", color = "Divergence from Ref")
dev.off()


# Fig5

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
                   aes(dist_to_tss, pval, color = rank), size = 2.25) +
        geom_point(data = best_df, 
                   aes(dist_to_tss, pval), 
                   shape = 1, size = 2.25, color = "black", stroke = 1) +
        geom_text_repel(data = best_df,
                        aes(dist_to_tss, pval, label = var_id),
                        fontface = "bold", size = 3, 
                        segment.color = "gray35",
                        nudge_x = 1.75e5, force = 1,
                        show.legend = FALSE) +
        scale_color_manual(values = c("0" = "#8491B4B2",
                                      "1" = "#DC0000B2",
                                      "2" = "gold3",
                                      "3" = "#7E6148FF",
                                      "4" = "#009E73")) +
        scale_x_continuous(breaks = seq(-5e5L, 5e5L, by = 2.5e5L),
                           labels = function(x) x/1e6L) +
        theme_bw() +
        theme(text = element_text(size = 10, family = "Arial"),
              strip.text = element_text(face = "bold"),
              strip.background = element_rect(fill = "grey")) +
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
    "./geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/2-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_mhc(mhc_genes)

ref_qtls <-
    "./geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/reference/3-conditional_analysis/conditional_60_all.txt.gz" %>%
    read_conditional_mhc(mhc_genes)

qtls_df <- 
    list("HLA-personalized" = hla_qtls, "Ref Transcriptome" = ref_qtls) %>%
    bind_rows(.id = "index") %>%
    mutate(gene = reorder(gene, tss))

hla_qtls <- hla_qtls %>%
    mutate(gene = reorder(gene, tss),
           colorize = ifelse(gene %in% gencode_hla$gene_name, 1L, 0L))

tiff("./plots_for_submission/Fig5.tiff",  width = 6.25, height = 7.25, units = "in", res = 300)
plot_qtls_1 <- ggplot() +
    geom_point(data = filter(hla_qtls, colorize == 0L), 
               aes(var_from, pval), 
               color = "grey", alpha = .1, size = .5) +
    geom_point(data = filter(hla_qtls, colorize == 1L),
               aes(var_from, pval, color = gene),
               alpha = .25, size = .5) +
    coord_cartesian(xlim = c(29.2e6L, 33.8e6L)) +
    ggsci::scale_color_npg() +
    scale_x_continuous(labels = function(x) x/1e6L) +
    theme(text = element_text(size = 9, family = "Arial"),
          axis.text = element_text(size = 8, family = "Arial"),
          axis.title.y = element_text(hjust = 0)) +
    guides(color = guide_legend(keyheight = .5, override.aes = list(alpha = 1, size = 3))) +
    labs(x = "position (Mb)",
         y = expression(paste("-log"[10], italic(Pvalue)))) 


plot_qtls_2 <- qtls_df %>%
    filter(gene %in% gencode_hla$gene_name) %>%
    plot_qtls_index()

plot_grid(NULL, plot_qtls_1, plot_qtls_2, ncol = 1,
          rel_heights = c(.025, .225, 1),
          labels = c("", "A", "B"), label_size = 12, vjust = 1, hjust = 0)
dev.off()


# Fig6


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
        scale_y_continuous(labels = comma, breaks = pretty_breaks(3)) +
        facet_wrap(~locus, scales = "free", ncol = 1, strip.position = "left") +
        labs(x = " ", y = "TPM") +
        theme_bw() +
        theme(text = element_text(size = 10, family = "Arial"),
              strip.text = element_text(face = "bold"),
              legend.position = "top") +
        guides(color = guide_legend(override.aes = list(size = 4)))
}

plot_slopes <- function(df) {
    ggplot(df, aes(reorder(id, resid, "mean"), resid)) +
        geom_jitter(width = .25, alpha = 1/2, size = .75) +
        geom_smooth(aes(group = 1), method = lm, se = FALSE) +
        scale_x_discrete(labels = function(x) sub("^([^_]+).+$", "\\1", x)) +
        scale_y_continuous(position = "right") +
        geom_text(data = distinct(df, locus, variant), 
                  aes(x = 1.5, y = 3, label = variant), size = 3) +
        coord_cartesian(ylim = c(-3, 3.2)) +
        facet_wrap(~locus, ncol = 1, scales = "free") +
        labs(x = " ", y = " ") +
        theme_bw() +
        theme(text = element_text(size = 10, family = "Arial"),
              strip.text = element_blank())
}

lineage_df <- 
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele, tpm) %>%
    mutate(subject = convert_ena_ids(subject),
           allele = gsub("IMGT_", "", allele),
           allele_3f = hla_trimnames(allele, 3),
           allele_2f = hla_trimnames(allele, 2),
           lineage = hla_trimnames(allele, 1))

dist_to_ref <- "./imgt_index/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)

haps_expression <- "./geuvadis_reanalysis/phase_hla/phase_hla_haps_snps.tsv" %>%
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
    "./geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/1-map_eqtls/th_50/1-phenotypes/phenotypes_60.bed.gz" %>%
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

best_eqtl_locus <- read_tsv("./geuvadis_reanalysis/eqtl_mapping/plots/eqtl.tsv") %>%
    filter(rank == 0) %>%
    select(locus = gene, variant = var_id)

eqtl_info <- "./geuvadis_reanalysis/phase_hla/eqtl_snps.vcf" %>%
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


tiff("./plots_for_submission/Fig6.tiff", width = 6, height = 8, units = "in", res = 300)
p1 <- plot_lineages(lineage_phased)

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- plot_slopes(eqtls_expression_df)

grid1 <- plot_grid(legend, NULL, p1, p2, ncol = 2, 
                   rel_widths = c(2.5, 1), rel_heights = c(.07, 1))

ggdraw(grid1) + 
    draw_label("HLA lineage", 0.45, 0.025, size = 10, fontfamily = "Arial") +
    draw_label("eQTL genotype", 0.82, 0.025, size = 10, fontfamily = "Arial") +
    draw_label("PCA-corrected expression", .985, 0.5, size = 10, angle = 90, fontfamily = "Arial")
dev.off()


# Fig7

calc_ase <- function(tpm) min(tpm)/sum(tpm)

ase_df <- 
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    filter(n_distinct(allele) == 2) %>%
    summarise(ase = calc_ase(tpm)) %>%
    ungroup()

tiff("./plots_for_submission/Fig7.tiff", width = 6, height = 3, units = "in", res = 300)
ggplot(ase_df, aes(reorder(locus, ase, median), ase)) +
    ggbeeswarm::geom_quasirandom(varwidth = TRUE, size = 1, alpha = 1/2) +
    stat_summary(fun.y = median, geom = "point", shape = "\U2014", size = 20, color = "#DC0000B2") +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(x = "", y = "ASE") +
    theme_bw() +
    theme(text = element_text(size = 10, family = "Arial"))
dev.off()


# Fig8

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

star_alleles <- 
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    mutate(subject = convert_ena_ids(subject),
           allele = gsub("IMGT_", "", allele)) %>%
    select(subject, locus, allele, tpm)

star_genes <- star_alleles %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

global_cors <- star_genes %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    spread(locus, tpm) %>% 
    select(-subject) %>% 
    ggcorr(label = TRUE, label_round = 2, label_size = 3.25, hjust = 0.6) +
    labs(title = "Gene-level") +
    theme(legend.position = "none",
          text = element_text(size = 9, family = "Arial"),
          plot.title = element_text(size = 11, family = "Arial"))


haps_data <- "./geuvadis_reanalysis/phase_hla/phase_hla_haps_snps.tsv" %>%
    read_tsv() %>% 
    distinct(subject, locus, hap, allele_gene) %>%
    rename(allele = allele_gene) %>%
    left_join(star_alleles, by = c("subject", "locus", "allele")) %>%
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
          plot.title = element_text(size = 11, family = "Arial"))


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
          plot.title = element_text(size = 11, family = "Arial"))


legend <- get_legend(trans_cors)

trans_cors <- trans_cors + theme(legend.position = "none")

plot_correlations <- 
    plot_grid(global_cors, cis_cors, trans_cors, legend, nrow = 1, 
              rel_widths = c(1, 1, 1, .25))


classII_and_CIITA <- gencode_chr_gene %>%
    filter(gene_name %in% c("HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "CIITA"))

class_2_trans_df <- 
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/gene_quantifications.tsv" %>%
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


tiff("./plots_for_submission/Fig8.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_grid(plot_correlations, plot_ciita, nrow = 2, ncol = 1, 
          rel_heights = c(1, .6), labels = c("A", "B"), label_size = 12, 
          hjust = 0)
dev.off()


# Fig S1

loci <- readLines("./imgt_index/imgt_loci.txt") %>%
    paste0("HLA-", .)

hla_personalized <- 
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% loci) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>% 
    left_join(geuvadis_info, by = c("subject" = "ena_id")) %>%
    select(subject = name, locus, tpm) %>%
    mutate(locus = reorder(locus, tpm, median),
           locus = factor(locus, levels = rev(levels(locus))),
           source = "HLA-personalized")

gencode12_hla <-
    "~/gencode_data/gencode.v12.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% loci)

phen <- 
    "./geuvadis_reanalysis/data/quantifications/peer/published/geuvadis_fpkms.csv" %>%
    read_csv()

hla_geuvadis <- phen %>%
    filter(subject %in% hla_personalized$subject) %>%
    select(subject, gencode12_hla$gene_id) %>%
    gather(gene_id, fpkm, -subject) %>%
    left_join(gencode12_hla, by = "gene_id") %>%
    select(subject, locus = gene_name, fpkm) %>%
    mutate(locus = factor(locus, levels = loci)) %>%
    complete(subject, locus, fill = list(fpkm = 0)) %>%
    mutate(locus = reorder(locus, fpkm, median),
           locus = factor(locus, levels = rev(levels(locus))),
           source = "Original Geuvadis")

p1 <- ggplot(hla_personalized, aes(locus, tpm)) +
    geom_boxplot() +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~source) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Arial"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = "TPM")

p2 <- ggplot(hla_geuvadis, aes(locus, fpkm)) +
    geom_boxplot() +
    geom_boxplot(fill = "grey") +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~source) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Arial"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = "FPKM")

tiff("./plots_for_submission/S1_fig.tiff", width = 6, height = 6, units = "in", res = 300)
plot_grid(p1, p2, ncol = 1)
dev.off()


# Fig S2

kallisto <- 
    "./geuvadis_reanalysis/expression/5-pseudoalignment/hla_personalized/quantifications_2/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup() %>%
    gather(unit, estimate, est_counts, tpm)

star_salmon <- 
    "./geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup() %>%
    gather(unit, estimate, est_counts, tpm)

star_kallisto_df <- 
    left_join(star_salmon, kallisto, by = c("subject", "locus", "unit"),
              suffix = c(".mapping", ".pseudo")) %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name))

cor_df_star_kallisto <- star_kallisto_df %>%
    group_by(locus, unit) %>%
    do(data.frame(r = cor(.$estimate.mapping, .$estimate.pseudo),
                  p = cor(.$estimate.mapping, .$estimate.pseudo, method = "spearman"),
                  x = min(.$estimate.mapping),
                  y = max(.$estimate.pseudo))) %>%
    ungroup() %>%
    mutate_at(vars(r, p), ~round(., digits = 2)) %>%
    mutate(label = paste("r == ", r, "*','~rho ==", p))

tiff("./plots_for_submission/S2_fig.tiff", width = 6, height = 8, units = "in", res = 300)
p_counts <- 
    ggplot(filter(star_kallisto_df, unit == "est_counts"),
           aes(estimate.mapping, estimate.pseudo)) +
    geom_abline() +
    geom_point(size = .8) +
    scale_x_continuous(breaks = pretty_breaks(n = 2), labels = comma) +
    scale_y_continuous(breaks = pretty_breaks(n = 3), labels = comma) +
    facet_wrap(~locus, scales = "free") +
    geom_text(data = filter(cor_df_star_kallisto, unit == "est_counts"), 
              aes(x, y, label = label),
              parse = TRUE, hjust = "inward", vjust = "inward", size = 3.5) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Arial"),
          axis.text = element_text(hjust = 1)) +
    labs(x = "STAR-Salmon", y = "kallisto", title = "Estimated Counts")

p_tpm <- 
    ggplot(filter(star_kallisto_df, unit == "tpm"),
           aes(estimate.mapping, estimate.pseudo)) +
    geom_abline() +
    geom_point(size = .8) +
    scale_x_continuous(breaks = pretty_breaks(n = 2), labels = comma) +
    scale_y_continuous(breaks = pretty_breaks(n = 3), labels = comma) +
    facet_wrap(~locus, scales = "free") +
    geom_text(data = filter(cor_df_star_kallisto, unit == "tpm"), 
              aes(x, y, label = label),
              parse = TRUE, hjust = "inward", vjust = "inward", size = 3.5) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Arial"),
          axis.text = element_text(hjust = 1)) +
    labs(x = "STAR-Salmon", y = "kallisto", title = "Transcripts per Million")

plot_grid(p_counts, NULL, p_tpm, nrow = 3, ncol = 1, rel_heights = c(1, 0.07, 1))
dev.off()


# Fig S3

gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_name, start, end, strand)

crd <- "./geuvadis_reanalysis/data/crd/LCL_ALL.chr6.subset.txt.gz" %>%
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

tiff("./plots_for_submission/S3_fig.tiff", width = 6, height = 3, units = "in", res = 300)
ggplot() +
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
    theme(text = element_text(size = 11, family = "Arial"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = c(.95, .8)) +
    labs(x = "Genomic position (Mb)", alpha = expression(~r^2)) +
    guides(alpha = guide_legend(override.aes = list(size = 2.5)))
dev.off()

# Fig S4

caveman <- 
    "./geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/caveman/results.hla" %>%
    read_tsv() %>%
    mutate(gene = factor(gene, levels = gencode_hla$gene_name),
           rank = as.character(rank))

tiff("./plots_for_submission/S4_fig.tiff", width = 7.5, height = 3, units = "in", res = 300)
ggplot(data = caveman, aes(rank, Probability, fill = index)) +
    geom_bar(stat = "identity", position = "dodge", alpha = .8) +
    ggsci::scale_fill_npg() +
    facet_wrap(~gene, nrow = 1, scales = "free_x") +
    labs(x = "") +
    theme(text = element_text(size = 11, family = "Arial"), 
          legend.position = "top")
dev.off()


# Fig S5
 
pcs <- seq(0, 100, 10)

egenes <- 
    "./geuvadis_reanalysis/eqtl_mapping/transcriptomemapping/hla_personalized/1-map_eqtls/th_50/2-permutations/results/permutations_%d.significant.txt" %>%
    sprintf(pcs) %>%
    setNames(pcs) %>%
    map_df(~read_delim(., delim = " ", col_names = FALSE), .id = "f") %>%
    count(f) %>%
    mutate(f = as.integer(f)) %>%
    arrange(f)

tiff("./plots_for_submission/S5_fig.tiff", width = 5, height = 3, units = "in", res = 300)
ggplot(egenes, aes(f, n)) + 
    geom_point(size = 2.5) + 
    geom_line() +
    scale_x_continuous(breaks = egenes$f) +
    scale_y_continuous(labels = comma) +
    scale_color_manual(values = cols) +
    theme_bw() +
    theme(text = element_text(size = 11, family = "Arial")) +
    labs(x = "Number of PCs", y = "Number of eGenes")
dev.off()


# Fig S6

make_pca_plot <- function(PC_x, PC_y) {
    
    ggplot(pcs, aes_string(PC_x, PC_y)) +
        geom_point(aes(color = pop)) +
        ggsci::scale_color_npg() +
        scale_x_continuous(breaks = pretty_breaks(3)) +
        scale_y_continuous(breaks = pretty_breaks(3)) +
        guides(color = guide_legend(override.aes = list(size = 2))) +
        theme(text = element_text(size = 11, family = "Arial")) +
        labs(color = "Population")
}

geuvadis_pops <- select(geuvadis_info, subject = name, pop)

pcs <- 
    read_delim("./geuvadis_reanalysis/data/pca_genotypes/eur.pca", delim = " ") %>%
    mutate(SampleID = sub("^.+_(PC\\d+)$", "\\1", SampleID)) %>%
    filter(SampleID %in% paste0("PC", 1:100)) %>%
    gather(subject, value, -1) %>%
    spread(SampleID, value) %>%
    mutate_at(vars(-subject), function(x) x/sqrt(sum(x^2))) %>%
    inner_join(geuvadis_pops, by = "subject")

tiff("./plots_for_submission/S6_fig.tiff", width = 7.5, height = 5, units = "in", res = 300)
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


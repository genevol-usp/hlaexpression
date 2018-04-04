devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(GGally)
library(scales)

# data
hla_genes <- gencode_hla$gene_name

geuvadis_ids <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(subject = ena_id, name)

# Boxplot
pipelines <- c("Ref Transcriptome", "HLA-personalized",
               "Ref Transcriptome (quasi)", "HLA-personalized (quasi)", 
               "Conventional", "Ref Genome (Unique)") 

cols <- ggsci::pal_npg()(6) %>%
    setNames(pipelines)

imgt_loci <- readLines("~/hlaexpression/imgt_index_v2/imgt_loci.txt") %>%
    paste0("HLA-", .)

mapping_imgt <- 
    "./3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% imgt_loci) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() %>%
    mutate(locus = reorder(locus, tpm, median),
           locus = factor(locus, levels = rev(levels(locus))))

mapping_ref <-
    "./3-map_to_transcriptome/reference/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    select(subject, locus, tpm)

quasi_imgt <- 
    "./4-quasimapping/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% imgt_loci) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup() 

quasi_ref <- 
    "./4-quasimapping/reference/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    select(subject, locus, tpm)

standard_imgt <- 
    read_tsv("./1-map_to_genome/quantifications/processed_imgt_quants.tsv") %>%
    select(subject, locus, tpm)

imgt_df <- list("Conventional" = standard_imgt, 
                "HLA-personalized" = mapping_imgt,
                "Ref Transcriptome" = mapping_ref,
                "HLA-personalized (quasi)" = quasi_imgt,
                "Ref Transcriptome (quasi)" = quasi_ref) %>% 
    bind_rows(.id = "pipeline") %>%
    group_by(locus) %>%
    filter(mean(tpm) >= 100) %>%
    ungroup() %>%
    mutate(locus = factor(locus, levels = levels(mapping_imgt$locus)),
           pipeline = factor(pipeline, levels = pipelines))

png("./plots/expression_boxplot.png", width = 8, height = 4, units = "in", res = 200)
ggplot(filter(imgt_df, pipeline == "HLA-personalized"), aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12)) +
    labs(y = "TPM")
dev.off()

png("./plots/expression_boxplot_pipelines.png", width = 10, height = 5, units = "in", res = 300)
ggplot(imgt_df, aes(locus, tpm)) +
    geom_boxplot(aes(fill = pipeline), outlier.size = .5) +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12),
          legend.position = "top") +
    labs(y = "TPM")
dev.off()

# TPM distributions
tpm_distribution_df <- imgt_df %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name))

png("./plots/tpm_distributions.png", height = 6, width = 10, units = "in", res = 200)
ggplot(tpm_distribution_df, aes(tpm, fill = pipeline)) +
    geom_density(alpha = .75) +
    scale_x_continuous(breaks = function(x) scales::pretty_breaks(3)(x)) +
    scale_fill_manual(values = cols) +
    facet_wrap(~locus, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "TPM")
dev.off()


# Scatterplots of pipeline comparisons
scatter_plot_cors <- function(df, x_var, y_var, dist_var) {
    
    cor_df <- df %>%
        group_by(locus) %>%
        do(data.frame(r = cor(.[[x_var]], .[[y_var]]),
                      p = cor(.[[x_var]], .[[y_var]], method = "spearman"),
                      x = min(.[[x_var]]),
                      y = max(.[[y_var]]))) %>%
        ungroup() %>%
        mutate_at(vars(r, p), ~round(., digits = 2)) %>%
        mutate(label = paste("r == ", r, "*','~rho ==", p))
        
    ggplot(df, aes_string(x_var, y_var)) +
        geom_abline() +
        geom_point(aes_string(color = dist_var), size = .8) +
        scale_color_gradient(low = "white", high = "darkred") +
        geom_text(data = cor_df, aes(x, y, label = label), parse = TRUE, 
                  hjust = "inward", vjust = "inward", size = 4, 
                  color = "white") +
        facet_wrap(~locus, scales = "free") +
        theme_dark() +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 16),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 16),
              legend.position = c(.75, .15))
}


## STAR vs kallisto
kallisto <- 
    "./5-pseudoalignment/hla_personalized/quantifications_2/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup() %>%
    gather(unit, estimate, est_counts, tpm)

star_salmon <- 
    "./3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
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

png("./plots/comparison_pseudoalignment.png", width = 6, height = 8, units = "in", res = 300)
p_counts <- 
    ggplot(filter(star_kallisto_df, unit == "est_counts"),
           aes(estimate.mapping, estimate.pseudo)) +
    geom_abline() +
    geom_point(size = .8) +
    scale_x_continuous(breaks = pretty_breaks(n = 2)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    facet_wrap(~locus, scales = "free") +
    geom_text(data = filter(cor_df_star_kallisto, unit == "est_counts"), 
              aes(x, y, label = label),
              parse = TRUE, hjust = "inward", vjust = "inward", size = 3.5) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, hjust = 1),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    labs(x = "STAR-Salmon", y = "kallisto",
         title = "Estimated Counts")

p_tpm <- 
    ggplot(filter(star_kallisto_df, unit == "tpm"),
           aes(estimate.mapping, estimate.pseudo)) +
    geom_abline() +
    geom_point(size = .8) +
    scale_x_continuous(breaks = pretty_breaks(n = 2)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    facet_wrap(~locus, scales = "free") +
    geom_text(data = filter(cor_df_star_kallisto, unit == "tpm"), 
              aes(x, y, label = label),
              parse = TRUE, hjust = "inward", vjust = "inward", size = 3.5) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, hjust = 1),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    labs(x = "STAR-Salmon", y = "kallisto", 
         title = "Transcripts per Million")

plot_grid(p_counts, NULL, p_tpm, nrow = 3, ncol = 1, rel_heights = c(1, 0.07, 1))
dev.off()


## Personalized vs reference
dist_ref <- read_tsv("../../imgt_index_v2/distances_to_reference.tsv")

mapping_imgt_dist <- 
    "./3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    left_join(dist_ref, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarise(dist = mean(dist)) %>%
    ungroup()

quasi_imgt_dist <- 
    "./4-quasimapping//hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    left_join(dist_ref, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarise(dist = mean(dist)) %>%
    ungroup()

dist_df <- left_join(mapping_imgt_dist, quasi_imgt_dist, 
                     by = c("subject", "locus"), suffix = c(".mapping", ".quasi"))

imgt_df_wide <- imgt_df %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    spread(pipeline, tpm) %>%
    rename(mapping.hla_personalized = `HLA-personalized`,
           mapping.ref = `Ref Transcriptome`,
           quasi.hla_personalized = `HLA-personalized (quasi)`,
           quasi.ref = `Ref Transcriptome (quasi)`,
           conventional = Conventional) %>%
    left_join(dist_df) %>%
    mutate(locus = factor(locus, levels = gencode_hla$gene_name))

### Mapping Pipeline
png("./plots/pers_vs_ref_mapping.png", height = 6, width = 10, units = "in", res = 300)
scatter_plot_cors(imgt_df_wide, "mapping.hla_personalized", "mapping.ref", "dist.mapping") +
    labs(x = "HLA-personalized", y = "Reference", color = "Divergence from Ref")
dev.off()


### Quasi-mapping Pipeline
png("./plots/pers_vs_ref_quasi.png", height = 6, width = 10, units = "in", res = 300)
scatter_plot_cors(imgt_df_wide, "quasi.hla_personalized", "quasi.ref", "dist.quasi") +
    labs(x = "HLA-personalized", y = "Reference", color = "Divergence from Ref")
dev.off()


# ASE
calc_ase <- function(tpm) min(tpm)/sum(tpm)

ase_df <- 
    "./3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    group_by(subject, locus) %>%
    filter(n_distinct(allele) == 2) %>%
    summarise(ase = calc_ase(tpm)) %>%
    ungroup()

png("./plots/ase.png", width = 8, height = 4, units = "in", res = 200)
ggplot(ase_df, aes(reorder(locus, ase, median), ase)) +
    ggbeeswarm::geom_quasirandom(varwidth = TRUE, size = 1, alpha = 1/2) +
    stat_summary(fun.y = median, geom = "point", shape = "\U2014", size = 20, color = "#DC0000B2") +
    coord_cartesian(ylim = c(0, 0.5)) +
    labs(x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))
dev.off()


# Correlation decrease as number of PCs increase
pcs <- c(seq(0, 20, 5), seq(30, 100, 10))

pca_star_df <-
    sprintf("../eqtl_mapping/transcriptomemapping/hla_personalized/1-phenotypes/phenotypes_eur_%d.bed.gz", pcs) %>%
    setNames(pcs) %>%
    map_df(. %>% read_tsv(progress = FALSE) %>% select(gid, matches("^HG|^NA")), 
           .id = "PC") %>%
    gather(subject, value, -(PC:gid))

hla_and_ciita <- filter(gencode_chr_gene, gene_name %in% c(hla_genes, "CIITA"))

pca_star_hla <- 
    inner_join(pca_star_df, hla_and_ciita, by = c("gid" = "gene_id")) %>%
    select(PC, subject, gene_name, value) %>%
    mutate(gene_name = sub("^HLA-", "", gene_name)) %>%
    spread(gene_name, value) 

cors_pca_star <- pca_star_hla %>%
    group_by(PC) %>%
    summarize(AxB = cor(A, B), 
              AxC = cor(A, C), 
              BxC = cor(B, C),
              DQA1xDQB1 = cor(DQA1, DQB1), 
              DQA1xDRB1 = cor(DQA1, DRB1), 
              DQB1xDRB1 = cor(DQB1, DRB1),  
              DQA1xCIITA = cor(DQA1, CIITA),
              DQB1xCIITA = cor(DQB1, CIITA), 
              DRB1xCIITA = cor(DRB1, CIITA)) %>%
    gather(gene_pair, correlation, -1) %>%
    mutate(PC = as.integer(PC)) %>%
    arrange(PC, gene_pair)

png("./plots/correlation_decrease.png", width = 10, height = 5, units = "in", res = 200)
ggplot(cors_pca_star, aes(PC, correlation)) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = pcs) +
    facet_wrap(~gene_pair) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.position = "top",	
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          axis.text.x = element_text(angle = 90)) +
    labs(x = "Number of PCs")
dev.off()

# Correlations
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
    "./3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
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
    ggcorr(label = TRUE, label_round = 2) +
    theme(legend.position = "none") +
    labs(title = "Gene-level")

haps_data <- "../phase_hla/phase_hla_haps_snps.tsv" %>%
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
    ggcorr(label = TRUE, label_round = 2) +
    theme(legend.position = "none") +
    labs(title = "Within haplotypes")

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
    ggcorr(data = NULL, cor_matrix = m, label = TRUE, label_round = 2, name = "r") +
    theme(legend.title = element_text(size = 14)) +
    labs(title = "Between haplotypes")

legend <- get_legend(trans_cors)

trans_cors <- trans_cors + theme(legend.position = "none")

png("./plots/correlations.png", width = 10, height = 4, units = "in", res = 300)
plot_grid(global_cors, cis_cors, trans_cors, legend, nrow = 1, 
          rel_widths = c(1, 1, 1, .2))
dev.off()

# Transactivator
classII_and_CIITA <- gencode_chr_gene %>%
    filter(gene_name %in% c("HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "CIITA"))

class_2_trans_df <- 
    "../expression/3-map_to_transcriptome/hla_personalized/quantifications_expressed50%.bed" %>%
    read_tsv() %>%
    inner_join(classII_and_CIITA, by = c("gid" = "gene_id")) %>%
    select(gid, locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, tpm, -gid, -locus) %>%
    select(subject, locus, tpm) %>%
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
    mutate(r = round(r, digits = 3))

png("./plots/trans_activ_corrs.png", width = 10, height = 3, units = "in", res = 200)
ggplot(class_2_trans_df, aes(tpm, CIITA)) +
    geom_point(size = 1) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(2)) +
    geom_text(data = cor_df, aes(x, y, label = paste("r =", r)),
              hjust = "inward", vjust = "inward", size = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 14, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16)) + 
    facet_wrap(~locus, nrow = 1, scales = "free") +
    labs(x = NULL)
dev.off()


# CRD
gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_name, start, end, strand)

crd <- "../../LCL_ALL.chr6.subset.txt.gz" %>%
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

png("./plots/crd.png", width = 6, height = 3, units = "in", res = 300)
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
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = c(.95, .8)) +
    labs(x = "Genomic position (Mb)", alpha = expression(~r^2)) +
    guides(alpha = guide_legend(override.aes = list(size = 2.5)))
dev.off()

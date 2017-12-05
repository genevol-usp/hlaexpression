devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(GGally)

# functions
calc_ase <- function(counts) min(counts)/sum(counts)

plot_lower <- function(data, mapping, ...) {
    
    ggplot(data = data, mapping = mapping) +
	geom_point(size = .5) +
	geom_smooth(method = lm) 
}

scatter_plot_cors <- function(df, x_var, y_var) {
    
    cor_df <- df %>%
        group_by(locus) %>%
        do(data.frame(r = cor(.[[x_var]], .[[y_var]]),
                      x = min(.[[x_var]]),
                      y = max(.[[y_var]]))) %>%
        ungroup() %>%
        mutate(r = round(r, digits = 3))
        
    ggplot(df, aes_string(x_var, y_var)) +
        geom_abline() +
        geom_point(size = .8) +
        geom_text(data = cor_df, aes(x, y, label = paste("r =", r)),
                  hjust = "inward", vjust = "inward", size = 5) +
        facet_wrap(~locus, scales = "free") +
        theme_bw() +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 16))
}

make_data <- function(df_allele, df_gene, locus1, locus2) {
  
    cis <- select(df_allele, subject, locus1, locus2)

    trans <- cis %>% 
	    group_by(subject) %>%
	    mutate_at(vars(!!locus2), rev) %>%
	    ungroup()

    gene <- select(df_gene, subject, locus1, locus2) %>%
        drop_na()

    bind_rows(list(Overall = gene, Cis = cis, Trans = trans), .id = "level") %>%
	    mutate(level = factor(level, levels = c("Overall", "Cis", "Trans")))
}

make_phase_plot <- function(data, locus1, locus2) {

    cor_df <- data %>%
        group_by(level) %>%
        do(data.frame(r = cor(.[[locus1]], .[[locus2]]),
                      x = min(.[[locus1]]),
                      y = max(.[[locus2]]))) %>%
        ungroup() %>%
        mutate(r = round(r, digits = 3))
    
    ggplot(data, aes_string(locus1, locus2)) +
	geom_point(size = .8) +
	geom_smooth(method = lm, se = FALSE) + 
    geom_text(data = cor_df, aes(x, y, label = paste("r =", r)),
              hjust = "inward", vjust = "inward", size = 5) +
	facet_wrap(~level, nrow = 1, scales = "free") +
	theme_bw() +
	theme(axis.text = element_text(size = 12),
	      axis.title = element_text(size = 16),
	      strip.text = element_text(size = 16)) + 
	labs(x = paste0("HLA-", locus1), y = paste0("HLA-", locus2))
}

# data
hla_genes <- sort(gencode_hla$gene_name)

geuvadis_ids <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(subject = ena_id, name)

## TPM
kallisto_imgt_tpm <- 
    read_tsv("./kallisto/quantifications_2/processed_quant.tsv") %>%
    filter(locus %in% hla_genes) %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, allele, est_counts, tpm) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup()

kallisto_pri_tpm <- 
    read_tsv("./kallisto/quantifications_PRI/processed_quant.tsv") %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, tpm) %>%
    arrange(subject, locus)

star_imgt <- 
    read_tsv("./star/quantifications_2/processed_quant.tsv") %>%
    filter(locus %in% hla_genes) %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, allele, est_counts, tpm) %>%
    mutate(allele = sub("IMGT_", "", allele))

star_imgt_tpm <-
    star_imgt %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm)) %>%
    ungroup()

star_imgt_with_nonclassical <- 
    read_tsv("./star/quantifications_2/processed_quant.tsv") %>%
    filter(grepl("HLA", locus)) %>%
    group_by(subject, locus) %>%
    summarize(tpm = sum(tpm)) %>%
    ungroup()

star_pri_tpm <- 
    read_tsv("./star/quantifications_PRI/processed_quant.tsv") %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, tpm)

## PCA
kallisto_imgt_pca <-
    "../qtls/qtls_kallisto/pca/imgt/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

kallisto_pri_pca <- 
    "../qtls/qtls_kallisto/pca/pri/phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_imgt_pca <-
    "../qtls/qtls_star/imgt/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_pri_pca <- 
    "../qtls/qtls_star/pri/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

quant_data <- 
    left_join(kallisto_imgt_tpm, star_imgt_tpm, by = c("subject", "locus"),
              suffix = c(".kallisto.imgt", ".star.imgt")) %>%
    left_join(kallisto_pri_tpm, by = c("subject", "locus")) %>%
    rename(tpm.kallisto.pri = tpm) %>%
    left_join(star_pri_tpm, by = c("subject", "locus")) %>%
    rename(tpm.star.pri = tpm) %>%
    left_join(kallisto_imgt_pca, by = c("subject", "locus")) %>%
    rename(pca.kallisto.imgt = resid) %>%
    left_join(star_imgt_pca, by = c("subject", "locus")) %>%
    rename(pca.star.imgt = resid) %>%
    left_join(kallisto_pri_pca, by = c("subject", "locus")) %>%
    rename(pca.kallisto.pri = resid) %>%
    left_join(star_pri_pca, by = c("subject", "locus")) %>%
    rename(pca.star.pri = resid)

star_tpm_df <- 
    quant_data %>%
    select(subject, locus, tpm.star.imgt, tpm.star.pri) %>%
    gather(index, tpm, 3:4) %>%
    mutate(index = sub("^tpm\\.star\\.", "", index)) %>%
    arrange(subject, locus, index)

pag3f <- mutate(pag, allele = hla_trimnames(allele, 3))

ase_df <- 
    star_imgt %>%
    group_by(subject, locus) %>%
    filter(n_distinct(allele) == 2) %>%
    summarize(ase = calc_ase(est_counts)) %>%
    ungroup()
  
ase_error <- 
    star_imgt %>%
    select(subject, locus, allele) %>%
    mutate(locus = sub("^HLA-", "", locus),
           allele = hla_trimnames(allele)) %>%
    calc_genotyping_accuracy(pag3f, by_locus = FALSE) %>%
    group_by(subject, locus) %>%
    summarize(error = sum(!correct)) %>%
    ungroup() %>%
    mutate(locus = paste0("HLA-", locus)) %>%
    inner_join(ase_df, by = c("subject", "locus"))

hla_and_transAct_genes <- gencode_chr_gene %>%
    filter(gene_name %in% c("HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "CIITA"))

class_2_trans_df <- 
    "../expression/star/quantifications_imgt_expressed90%.csv" %>%
    read_csv() %>%
    select(subject, hla_and_transAct_genes$gene_id) %>%  
    gather(gene, tpm, -subject) %>%
    left_join(hla_and_transAct_genes, by = c("gene" = "gene_id")) %>%
    select(subject, locus = gene_name, tpm) %>%
    spread(locus, tpm) %>%
    gather(locus, tpm, -subject, -CIITA) %>%
    select(subject, locus, tpm, CIITA) %>%
    arrange(subject, locus)

pcs <- c(seq(0, 20, 5), seq(30, 100, 10))

pca_star_df <-
    sprintf("../qtls/qtls_star/imgt/1-phenotypes/phenotypes_eur_%d.bed.gz", pcs) %>%
    setNames(pcs) %>%
    map_df(. %>% read_tsv(progress = FALSE) %>% select(gid, matches("^HG|^NA")), 
	   .id = "PC") %>%
    gather(subject, value, -(PC:gid))

pca_star_hla <-
    inner_join(pca_star_df, hla_and_transAct_genes, by = c("gid" = "gene_id")) %>%
    select(PC, subject, gene_name, value) %>%
    mutate(gene_name = sub("^HLA-", "", gene_name)) %>%
    spread(gene_name, value) 

cors_pca_star <- 
    pca_star_hla %>%
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

concordant_haps_I <- 
    read_tsv("./star/phase_hla_alleles/data/concordant_haps_classI.tsv") %>%
    select(subject, hap)

concordant_haps_II <- 
    read_tsv("./star/phase_hla_alleles/data/concordant_haps_classII.tsv") %>%
    select(subject, hap)

tpm_by_allele_wide_I <- 
    read_tsv("./star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv") %>%
    select(subject, locus, hap, tpm) %>%
    filter(locus %in% c("A", "B", "C")) %>%
    spread(locus, tpm) %>%
    inner_join(concordant_haps_I, by = c("subject", "hap")) %>%
    arrange(subject, hap)

tpm_by_allele_wide_II <- 
    read_tsv("./star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv") %>%
    select(subject, locus, hap, tpm) %>%
    filter(!locus %in% c("A", "B", "C")) %>%
    spread(locus, tpm) %>%
    inner_join(concordant_haps_II, by = c("subject", "hap")) %>%
    arrange(subject, hap)

tpm_by_gene_wide <- 
    full_join(tpm_by_allele_wide_I, tpm_by_allele_wide_II, by = c("subject", "hap")) %>%
    group_by(subject) %>%
    summarize_at(vars(A:DRB1), sum) %>%
    ungroup()

# plots
png("./plots/star_vs_kallisto_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "tpm.star.imgt", "tpm.kallisto.imgt") +
    labs(x = "TPM (STAR-Salmon)", y = "TPM (kallisto)")
dev.off()

png("./plots/star_vs_kallisto_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "resid.star.imgt.tpm", "resid.kallisto.imgt.tpm") +
    labs(x = "PCA-corrected estimate (STAR-Salmon)", y = "PCA-corrected estimate (kallisto)")
dev.off()

png("./plots/star_imgt_vs_pri_TPM.png", height = 6, width = 10, units = "in", res = 200)
scatter_plot_cors(quant_data, "tpm.star.imgt", "tpm.star.pri") +
    labs(x = "TPM (HLA personalized index)", y = "TPM (Ref transcriptome)")
dev.off()

png("./plots/star_imgt_vs_pri_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "resid.star.imgt.tpm", "resid.star.pri.tpm") +
    labs(x = "PCA-corrected estimate (HLA personalized index)", 
         y = "PCA-corrected estimate (Ref transcriptome)")
dev.off()

png("./plots/kallisto_imgt_vs_pri_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "tpm.kallisto.imgt", "tpm.kallisto.pri") +
    labs(x = "TPM (HLA personalized index)", y = "TPM (Ref transcriptome)")
dev.off()

png("./plots/kallisto_imgt_vs_pri_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_cors(quant_data, "resid.kallisto.imgt.tpm", "resid.kallisto.pri.tpm") +
    labs(x = "PCA-corrected estimate (HLA personalized index)", 
         y = "PCA-corrected estimate (Ref transcriptome)")
dev.off()

png("./plots/tpm_distributions.png", height = 6, width = 10, units = "in", res = 200)
ggplot(star_tpm_df, aes(tpm, fill = index)) +
    geom_density(alpha = 1/2) +
    scale_x_continuous(breaks = function(x) scales::pretty_breaks(3)(x)) +
    ggthemes::scale_fill_colorblind() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, margin = margin(t = 10))) +
    facet_wrap(~locus, scales = "free")
dev.off()

png("./plots/ase.png", width = 8, height = 5, units = "in", res = 200)
ggplot(ase_error, aes(factor(error), ase)) +
    ggbeeswarm::geom_quasirandom(varwidth = TRUE, size = .75, alpha = 1/2) +
    scale_y_continuous(limits = c(0, 0.5)) +
    facet_wrap(~locus) + 
    labs(x = "number of wrong calls in genotype") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
	  axis.title = element_text(size = 16),
	  strip.text = element_text(size = 16))
dev.off()

png("./plots/ase_histogram.png", width = 8, height = 4, units = "in", res = 200)
ggplot(ase_df, aes(ase)) +
  geom_density(fill = "grey35", color = NA) +
  facet_wrap(~locus) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))
dev.off()

png("./plots/hlacorrelations.png", width = 6, height = 6, units = "in", res = 200)
pairs_hla_k <-
    star_imgt_tpm %>%
    select(-est_counts) %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    spread(locus, tpm) %>%
    select(-subject) %>%
    ggpairs(lower = list(continuous = plot_lower), upper = list()) + 
    theme_bw() +
    theme(title = element_text(size = 14),
          axis.text.x = element_text(angle = 90))

print(pairs_hla_k, left = .3, bottom = .3)
dev.off()

png("./plots/trans_activ_corrs.png", width = 10, height = 3.5, units = "in", res = 200)
cor_df <- class_2_trans_df %>%
    group_by(locus) %>%
    do(data.frame(r = cor(.$tpm, .$CIITA),
                  x = min(.$tpm),
                  y = max(.$CIITA))) %>%
    ungroup() %>%
    mutate(r = round(r, digits = 3))

ggplot(class_2_trans_df, aes(tpm, CIITA)) +
    geom_point(size = .8) +
    geom_smooth(method = lm, se = FALSE) +
    scale_x_continuous(breaks = scales::pretty_breaks(2)) +
    geom_text(data = cor_df, aes(x, y, label = paste("r =", r)),
              hjust = "inward", vjust = "inward", size = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 14, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16)) + 
    facet_wrap(~locus, nrow = 1) +
    labs(x = NULL)
dev.off()

png("./plots/a_vs_b.png", height = 3.5, width = 10, units = "in", res = 200)
a_b <- make_data(tpm_by_allele_wide_I, tpm_by_gene_wide, "A", "B")
make_phase_plot(a_b, "A", "B")
dev.off()

png("./plots/a_vs_c.png", height = 3.5, width = 10, units = "in", res = 200)
a_c <- make_data(tpm_by_allele_wide_I, tpm_by_gene_wide, "A", "C")
make_phase_plot(a_c, "A", "C")
dev.off()

png("./plots/b_vs_c.png", height = 3.5, width = 10, units = "in", res = 200)
b_c <- make_data(tpm_by_allele_wide_I, tpm_by_gene_wide, "B", "C")
make_phase_plot(b_c, "B", "C")
dev.off()

png("./plots/dqa_vs_dqb.png", height = 3.5, width = 10, units = "in", res = 200)
dqa_dqb <- make_data(tpm_by_allele_wide_II, tpm_by_gene_wide, "DQA1", "DQB1")
make_phase_plot(dqa_dqb, "DQA1", "DQB1")
dev.off()

png("./plots/dqa_vs_drb.png", height = 3.5, width = 10, units = "in", res = 200)
dqa_drb <- make_data(tpm_by_allele_wide_II, tpm_by_gene_wide, "DQA1", "DRB1")
make_phase_plot(dqa_drb, "DQA1", "DRB1")
dev.off()

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

png("./plots/expression_boxplot.png", width = 8, height = 5, units = "in", res = 200)
ggplot(star_imgt_with_nonclassical, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 12)) +
    labs(y = "TPM")
dev.off()

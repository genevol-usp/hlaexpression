devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
library(cowplot)
library(GGally)
library(ggpmisc)

# functions
calc_ase <- function(counts) min(counts)/sum(counts)

plot_lower <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
  geom_point(size = .6, alpha = .5) +
  geom_smooth(method = lm) 
}

scatter_plot <- function(df, x_var, y_var) {
  ggplot(df, aes_string(x_var, y_var)) +
    geom_abline() +
    geom_point() +
    facet_wrap(~locus, scales = "free") +
    ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
                          formula = y ~ x, parse = TRUE, size = 6) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))
}

make_data <- function(locus1, locus2) {
  l1 <- enquo(locus1)
  l2 <- enquo(locus2)

  cis <- select(tpm_by_allele_wide, subject, !!l1, !!l2)
  trans <- cis %>% 
    group_by(subject) %>%
    mutate(!!quo_name(l2) := rev(!!l2)) %>%
    ungroup()
  gene <- select(tpm_by_gene_wide, subject, !!l1, !!l2)
  
  bind_rows(list(Overall = gene, Cis = cis, Trans = trans), .id = "level") %>%
    mutate(level = factor(level, levels = c("Overall", "Cis", "Trans")))
}

make_phase_plot <- function(data, locus1, locus2) {
  ggplot(data, aes_string(locus1, locus2)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) + 
  facet_wrap(~level, nrow = 1, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) + 
  labs(x = paste0("HLA-", locus1), y = paste0("HLA-", locus2)) +
  stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
               formula = y ~ x, parse = TRUE, size = 6)
}

# data
hla_genes <- paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1")) 

geuvadis_ids <- geuvadis_info %>%
  filter(kgp_phase3 == 1, pop != "YRI") %>%
  select(subject = ena_id, name)

gencode_hla <- gencode_chr_gene %>%
  filter(gene_name %in% hla_genes)

distances <- read_tsv("../../simulation/data/distances_to_reference.tsv")

kallisto_imgt_tpm <- 
    read_tsv("./kallisto/quantifications_2/processed_quant.tsv") %>%
    filter(locus %in% hla_genes) %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, allele, est_counts, tpm) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    left_join(distances, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm), dist = mean(dist)) %>%
    ungroup()

kallisto_pri_tpm <- 
    read_tsv("./kallisto/quantifications_PRI/processed_quant.tsv") %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, tpm) %>%
    arrange(subject, locus)

kallisto_imgt_10pc <-
    "../qtls/qtls_kallisto/pca_correction/imgt/phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

kallisto_pri_10pc <- 
    "../qtls/qtls_kallisto/pca_correction/pri/phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_imgt <- 
    read_tsv("./star/quantifications_2/processed_quant.tsv") %>%
    filter(locus %in% hla_genes) %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, allele, est_counts, tpm) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    left_join(distances, by = c("locus", "allele")) 

star_imgt_tpm <-
    star_imgt %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm), dist = mean(dist)) %>%
    ungroup()

star_pri_tpm <- 
    read_tsv("./star/quantifications_PRI/processed_quant.tsv") %>%
    inner_join(geuvadis_ids, by = "subject") %>%
    select(subject = name, locus, tpm)

star_imgt_10pc <-
    "../qtls/qtls_star/imgt/1-phenotypes/phenotypes_eur_10.bed.gz" %>%
    read_tsv() %>%
    inner_join(select(gencode_hla, gene_id, gene_name), by = c("gid" = "gene_id")) %>%
    select(locus = gene_name, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, resid, -locus)

star_pri_10pc <- 
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
    left_join(kallisto_imgt_10pc, by = c("subject", "locus")) %>%
    rename(resid.kallisto.imgt.tpm = resid) %>%
    left_join(star_imgt_10pc, by = c("subject", "locus")) %>%
    rename(resid.star.imgt.tpm = resid) %>%
    left_join(kallisto_pri_10pc, by = c("subject", "locus")) %>%
    rename(resid.kallisto.pri.tpm = resid) %>%
    left_join(star_pri_10pc, by = c("subject", "locus")) %>%
    rename(resid.star.pri.tpm = resid)

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

hla_and_transAct_genes <- 
    filter(gencode_chr_gene, gene_name %in% c(hla_genes, "CIITA"))

pcs <- c(seq(0, 20, 5), seq(30, 100, 10))

pca_star_df <-
    sprintf("../qtls/qtls_star/imgt/1-phenotypes/phenotypes_eur_%d.bed.gz", pcs) %>%
    setNames(pcs) %>%
    parallel::mclapply(function(i) 
        read_tsv(i, progress = FALSE) %>% 
            inner_join(hla_and_transAct_genes, by = c("id" = "gene_id")) %>%
            select(gene_name, HG00096:NA20828),
        mc.cores = length(pcs)) %>%
    bind_rows(.id = "PCs") %>%
    gather(subject, value, HG00096:NA20828) %>%
    mutate(gene_name = sub("^HLA-", "", gene_name)) %>%
    spread(gene_name, value) 

phen10 <- filter(pca_star_df, PCs == 10) %>% select(-PCs)

class_2_trans_df <- phen10 %>%
    select(subject, DQA1, DQB1, DRB1, CIITA) %>%
    gather(locus, value, DQA1, DQB1, DRB1) %>%
    arrange(subject, locus)

concordant_haps <- 
    read_tsv("./star/phase_hla_alleles/data/concordant_haps_classIandII.tsv") %>%
    select(subject, hap)

tpm_by_allele_wide <- 
    read_tsv("./star/phase_hla_alleles/data/1000G_haps_expression_snps.tsv") %>%
    filter(rank == 0L) %>%
    select(subject, locus, hap, tpm) %>%
    spread(locus, tpm) %>%
    inner_join(concordant_haps, by = c("subject", "hap")) %>%
    arrange(subject, hap)

tpm_by_gene_wide <- 
    star_imgt_tpm %>%
    mutate(locus = sub("HLA-", "", locus)) %>%
    select(subject, locus, tpm) %>%
    spread(locus, tpm) %>%
    filter(subject %in% concordant_haps$subject)

cors_pca_star <- 
    pca_star_df %>%
    group_by(PCs) %>%
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
    mutate(PCs = as.integer(PCs)) %>%
    arrange(PCs, gene_pair)

# plots
png("./plots/star_vs_kallisto_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.kallisto.imgt", "tpm.star.imgt") +
    labs(x = "TPM (kallisto)", y = "TPM (STAR-Salmon)")
dev.off()

png("./plots/star_vs_kallisto_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "resid.kallisto.imgt.tpm", "resid.star.imgt.tpm") +
    labs(x = "PCA-corrected TPM (kallisto)", y = "PCA-corrected TPM (STAR-Salmon)")
dev.off()

png("./plots/star_imgt_vs_pri_TPM.png", height = 6, width = 10, units = "in", res = 200)
scatter_plot(quant_data, "tpm.star.imgt", "tpm.star.pri") +
    labs(x = "TPM (STAR-IMGT)", y = "TPM (STAR REF chromosomes)")
dev.off()

png("./plots/star_imgt_vs_pri_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "resid.star.imgt.tpm", "resid.star.pri.tpm") +
    labs(x = "PCA-corrected TPM (STAR-IMGT)", 
         y = "PCA-corrected TPM (STAR REF chromosomes)")
dev.off()

png("./plots/kallisto_imgt_vs_pri_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.kallisto.imgt", "tpm.kallisto.pri") +
    labs(x = "TPM (kallisto-IMGT)", y = "TPM (kallisto REF chromosomes)")
dev.off()

png("./plots/kallisto_imgt_vs_pri_PCA.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "resid.kallisto.imgt.tpm", "resid.kallisto.pri.tpm") +
    labs(x = "PCA-corrected TPM (kallisto-IMGT)", 
         y = "PCA-corrected TPM (kallisto REF chromosomes)")
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

png("./plots/hlacorrelations.png", width = 8, height = 8, units = "in", res = 200)
pairs_hla_k <-
  ggpairs(select(phen10, -subject, -CIITA), 
          lower = list(continuous = plot_lower), upper = list()) + 
  theme_bw() +
  theme(title = element_text(size = 14))

print(pairs_hla_k, left = .3, bottom = .3)
dev.off()

png("./plots/trans_activ_corrs.png", width = 10, height = 3.5, units = "in", res = 200)
ggplot(class_2_trans_df, aes(value, CIITA)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16)) + 
  facet_wrap(~locus) +
  labs(x = NULL) +
  stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
               formula = y ~ x, parse = TRUE, size = 6)
dev.off()

png("./plots/a_vs_b.png", height = 3.5, width = 10, units = "in", res = 200)
a_b <- make_data(A, B)
make_phase_plot(a_b, "A", "B")
dev.off()

png("./plots/a_vs_c.png", height = 3.5, width = 10, units = "in", res = 200)
a_c <- make_data(A, C)
make_phase_plot(a_c, "A", "C")
dev.off()

png("./plots/b_vs_c.png", height = 3.5, width = 10, units = "in", res = 200)
b_c <- make_data(B, C)
make_phase_plot(b_c, "B", "C")
dev.off()

png("./plots/dqa_vs_dqb.png", height = 3.5, width = 10, units = "in", res = 200)
dqa_dqb <- make_data(DQA1, DQB1)
make_phase_plot(dqa_dqb, "DQA1", "DQB1")
dev.off()

png("./plots/dqa_vs_drb.png", height = 3.5, width = 10, units = "in", res = 200)
dqa_drb <- make_data(DQA1, DRB1)
make_phase_plot(dqa_drb, "DQA1", "DRB1")
dev.off()

png("./plots/correlation_decrease.png", width = 10, height = 5, units = "in", res = 200)
ggplot(cors_pca_star, aes(PCs, correlation)) +
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
  labs(x = "Number of PCs/factors")
dev.off()

png("./plots/expression_boxplot.png", width = 8, height = 5, units = "in", res = 200)
ggplot(star_imgt_tpm, aes(locus, tpm)) +
    geom_boxplot(fill = "grey") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.text = element_text(size = 12)) +
    labs(y = "TPM")
dev.off()

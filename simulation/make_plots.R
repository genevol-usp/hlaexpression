devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

# Functions
read_pcanorm <- function(f) {

    read_tsv() %>%
	inner_join(select(gencode_hla, gene_id, gene_name), 
		   by = c("gid" = "gene_id")) %>%
	select(locus = gene_name, starts_with("sample_")) %>%
	gather(subject, resid, -locus) %>%
	mutate(locus = sub("^HLA-", "", locus))
}

plot_dist <- function(df) {

    ggplot(df, aes(dist, prop_mapped, color = index)) +
	geom_point() +
	geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
		  alpha = 0.4, size = 1.5) +
	scale_x_continuous(labels = scales::percent) +
	scale_y_continuous(breaks = seq(0, 1.5, 0.5)) +
	facet_wrap(~locus, scales = "free_x") +
	ggsci::scale_color_aaas(labels = c(pri = "Ref chromosomes",
					  imgt = stringr::str_wrap("Ref chromosomes + HLA diversity", 20))) +
	theme_bw() +
	theme(axis.text = element_text(size = 14),
	      axis.title = element_text(size = 18),
	      legend.title = element_text(size = 16),
	      legend.text = element_text(size = 14),
	      strip.text = element_text(size = 16, face = "bold"),
	      legend.position = "top") +
	guides(color = guide_legend(override.aes = list(size = 4))) +
	labs(x = "proportion of sites with mismatches to the reference allele", 
	     y = "proportion of reads recovered")
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

scatter_plot_color <- function(df, x_var, y_var, alpha_var) {

    ggplot(df, aes_string(x_var, y_var, alpha = alpha_var)) +
	geom_abline() +
	geom_point() +
	facet_wrap(~locus, scales = "free") +
	ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 2,
			      formula = y ~ x, parse = TRUE, size = 6) +
	theme_bw() +
	theme(axis.text = element_text(size = 12),
	      axis.title = element_text(size = 16),
	      strip.text = element_text(size = 16)) +
	labs(alpha = "distance to ref")
}

# Data
allele_dist <- read_tsv("./data/distances_to_reference.tsv") %>%
    mutate(locus = sub("^HLA-", "", locus))

samples <- tibble(subject = readLines("./data/samples.txt"),
                  code = sprintf("sample_%02d", 1:50))

gencode_hla <- gencode_chr_gene %>%
    filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1")))

index <- Biostrings::readDNAStringSet("./data/polyester_index.fa")

ground_truth <- 
    read_tsv("./data/phenotypes.tsv") %>%
    mutate(target_id = names(index)) %>%
    filter(grepl("IMGT_(A|B|C|DQA1|DQB1|DRB1)", target_id)) %>%
    gather(subject, true_counts, -target_id) %>%
    inner_join(samples, by = "subject") %>%
    select(subject = code, target_id, true_counts) %>%
    mutate(locus = sub("^IMGT_([^*]+).+$", "\\1", target_id)) %>%
    group_by(subject, locus) %>%
    summarize(true_counts = sum(true_counts))

kallisto_quant_imgt <- 
    read_tsv("./expression/kallisto/quantifications_2/processed_quant.tsv") %>%
    filter(locus %in% c("A", "B", "C", "DQA1", "DQB1", "DRB1")) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    left_join(allele_dist, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm), dist = mean(dist)) %>%
    ungroup()

kallisto_quant_pri <- 
    read_tsv("./expression/kallisto/quantifications_PRI/processed_quant.tsv")

star_quant_imgt <- 
    read_tsv("./expression/star/quantifications_2/processed_quant.tsv") %>%
    filter(locus %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))) %>%
    mutate(locus = sub("HLA-", "", locus), allele = sub("IMGT_", "", allele)) %>%
    left_join(allele_dist, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), tpm = sum(tpm), dist = mean(dist)) %>%
    ungroup()

star_quant_pri <- 
    read_tsv("./expression/star/quantifications_PRI/processed_quant.tsv") %>%
    mutate(locus = sub("HLA-", "", locus))

kallisto_10pc <-
    "./expression/kallisto/phenotype_correction/imgt/phenotypes_10.bed.gz" %>%
    read_pcanorm()

star_10pc <-
    "./expression/star/phenotype_correction/imgt/phenotypes_10.bed.gz" %>%
    read_pcanorm()

kallisto_pri_10pc <-
    "./expression/kallisto/phenotype_correction/pri/phenotypes_10.bed.gz" %>%
    read_pcanorm()

star_pri_10pc <-
    "./expression/star/phenotype_correction/priphenotypes_10.bed.gz" %>%
    read_pcanorm()

quant_data <-
    left_join(kallisto_quant_imgt, kallisto_quant_pri, 
	      by = c("subject", "locus"), 
	      suffix = c(".kallisto.imgt", ".kallisto.pri")) %>%
    left_join(star_quant_imgt, by = c("subject", "locus")) %>%
    rename(est_counts.star.imgt = est_counts, tpm.star.imgt = tpm,
	   dist.star = dist) %>%
    left_join(star_quant_pri, by = c("subject", "locus")) %>%
    rename(est_counts.star.pri = est_counts, tpm.star.pri = tpm) %>%
    left_join(kallisto_10pc, by = c("subject", "locus")) %>%
    rename(resid.kallisto.imgt = resid) %>%
    left_join(star_10pc, by = c("subject", "locus")) %>%
    rename(resid.star.imgt = resid) %>%
    left_join(kallisto_pri_10pc, by = c("subject", "locus")) %>%
    rename(resid.kallisto.pri = resid) %>%
    left_join(star_pri_10pc, by = c("subject", "locus")) %>%
    rename(resid.star.pri = resid)

# Plots
png("./plots/kallisto_prop_mapped.png", width = 12, height = 6, units = "in", res = 200)
plot_dist(counts_kallisto)
dev.off()

png("./plots/star_prop_mapped.png", width = 8, height = 5, units = "in", res = 200)
plot_dist(counts_star)
dev.off()

png("./plots/kallisto_vs_star_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.kallisto.imgt", "tpm.star.imgt") +
  labs(x = "TPM (kallisto)", 
       y = "TPM (STAR-Salmon)")
dev.off()

png("./plots/kallisto_vs_star_10pc.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "resid.kallisto.imgt", "resid.star.imgt") +
  labs(x = "PCA-corrected TPM (kallisto)", 
       y = "PCA-corrected TPM (STAR-Salmon)")
dev.off()

png("./plots/kallisto_vs_star_PRI_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.kallisto.pri", "tpm.star.pri") +
  labs(x = "TPM (kallisto)", 
       y = "TPM (STAR-Salmon)")
dev.off()

png("./plots/kallisto_vs_star_PRI_10pc.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "resid.kallisto.pri", "resid.star.pri") +
  labs(x = "PCA-corrected TPM (kallisto)", 
       y = "PCA-corrected TPM (STAR-Salmon)")
dev.off()

png("./plots/kallisto_imgt_vs_chr_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_color(quant_data, "tpm.kallisto.imgt", "tpm.kallisto.pri", "dist.kallisto") +
  labs(x = "TPM (kallisto-IMGT)", 
       y = "TPM (kallisto REF chromosomes)")
dev.off()

png("./plots/kallisto_imgt_vs_pri_10pc.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_color(quant_data, "resid.kallisto.imgt", "resid.kallisto.pri", "dist.kallisto") +
  labs(x = "PCA-corrected TPM (kallisto-IMGT)", 
       y = "PCA-corrected TPM (kallisto REF chromosomes)")
dev.off()

png("./plots/star_imgt_vs_pri_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_color(quant_data, "tpm.star.imgt", "tpm.star.pri", "dist.star") +
  labs(x = "TPM (STAR-IMGT)", 
       y = "TPM (STAR REF chromosomes)")
dev.off()

png("./plots/star_imgt_vs_chr_10pc.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot_color(quant_data, "resid.star.imgt", "resid.star.pri", "dist.star") +
  labs(x = "PCA-corrected TPM (STAR-IMGT)", 
       y = "PCA-corrected TPM (STAR REF chromosomes)")
dev.off()

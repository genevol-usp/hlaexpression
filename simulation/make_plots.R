devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

# Functions
read_imgt_quants <- function(f) {
  
    read_tsv(f) %>%
        filter(locus %in% hla_genes) %>%
        mutate(allele = gsub("IMGT_", "", allele)) %>%
        left_join(allele_dist, by = c("locus", "allele")) %>%
        group_by(subject, locus) %>%
        summarize(est_counts = sum(est_counts), tpm = sum(tpm), 
                  dist = mean(dist)) %>%
        ungroup()
}
 
read_pcanorm <- function(f) {

    read_tsv(f) %>%
	    inner_join(select(gencode_hla, gene_id, gene_name), 
	               by = c("gid" = "gene_id")) %>%
	    select(locus = gene_name, starts_with("sample_")) %>%
	    gather(subject, resid, -locus)
}

plot_dist <- function(df) {

    ggplot(df, aes(dist, prop_mapped, color = index)) +
	    geom_point() +
	    geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
	              alpha = 0.4, size = 1.5) +
	    scale_x_continuous(labels = scales::percent,
	                       breaks = scales::pretty_breaks(n = 3)) +
    	scale_y_continuous(breaks = seq(0, 1.5, 0.5)) +
    	facet_wrap(~locus, scales = "free_x") +
    	scale_color_manual(
    	    values = c(imgt = "#8491B4B2", pri = "#DC0000B2"),
    	    labels = c(pri = "Ref transcriptome",
    	               imgt = "Personalized index")) +
    	theme_bw() +
    	theme(axis.text = element_text(size = 8),
    	      axis.title = element_text(size = 10),
    	      legend.title = element_text(size = 12),
    	      legend.text = element_text(size = 10),
    	      strip.text = element_text(size = 10, face = "bold"),
    	      legend.position = c(.75, .15)) +
    	guides(color = guide_legend(override.aes = list(size = 4))) +
    	labs(x = "sequence divergence to the HLA reference allele", 
    	     y = "estimated counts / ground truth")
}

scatter_plot <- function(df, x_var, y_var) {
  
    ggplot(df, aes_string(x_var, y_var)) +
	geom_abline() +
    	geom_point() +
    	facet_wrap(~locus, scales = "free") +
    	ggpmisc::stat_poly_eq(aes(label = ..adj.rr.label..), rr.digits = 3,
    	                      formula = y ~ x, parse = TRUE, size = 6) +
    	theme_bw() +
    	theme(axis.text = element_text(size = 12),
    	      axis.title = element_text(size = 16),
    	      strip.text = element_text(size = 16))
}

# Data
allele_dist <- read_tsv("./PEreads_75bp/data/distances_to_reference.tsv")

hla_genes <- sort(gencode_hla$gene_name)

index <- Biostrings::readDNAStringSet("./PEreads_75bp/data/polyester_index.fa")

ground_truth <- read_tsv("./PEreads_75bp/data/phenotypes.tsv") %>%
    mutate(target_id = names(index)) %>%
    filter(grepl("IMGT_(A|B|C|DPB1|DQA1|DQB1|DRB1)", target_id)) %>%
    gather(subject, true_counts, -target_id) %>%
    mutate(locus = sub("^IMGT_([^*]+).+$", "HLA-\\1", target_id)) %>%
    group_by(subject, locus) %>%
    summarize(true_counts = sum(true_counts))

kallisto_quant_imgt <- 
    "./PEreads_75bp/expression/kallisto/quantifications_2/processed_quant.tsv" %>%
    read_imgt_quants()

kallisto_quant_pri <- 
    read_tsv("./PEreads_75bp/expression/kallisto/quantifications_PRI/processed_quant.tsv")

star_quant_imgt <- 
    "./PEreads_75bp/expression/star/quantifications_2/processed_quant.tsv" %>%
    read_imgt_quants()

star_quant_pri <- 
    read_tsv("./PEreads_75bp/expression/star/quantifications_PRI/processed_quant.tsv")

kallisto_10pc <-
    "./PEreads_75bp/expression/kallisto/phenotype_correction/imgt/phenotypes_10.bed.gz" %>%
    read_pcanorm()

star_10pc <-
    "./PEreads_75bp/expression/star/phenotype_correction/imgt/phenotypes_10.bed.gz" %>%
    read_pcanorm()

kallisto_pri_10pc <-
    "./PEreads_75bp/expression/kallisto/phenotype_correction/pri/phenotypes_10.bed.gz" %>%
    read_pcanorm()

star_pri_10pc <-
    "./PEreads_75bp/expression/star/phenotype_correction/pri/phenotypes_10.bed.gz" %>%
    read_pcanorm()

quant_data <-
    left_join(kallisto_quant_imgt, star_quant_imgt, 
	      by = c("subject", "locus"), 
	      suffix = c(".kallisto.imgt", ".star.imgt")) %>%
    left_join(kallisto_quant_pri, by = c("subject", "locus")) %>%
    rename(est_counts.kallisto.pri = est_counts, tpm.kallisto.pri = tpm) %>%
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

counts_kallisto <- quant_data %>%
    select(subject, locus, dist.kallisto.imgt, starts_with("est_counts.kallisto")) %>%
    gather(index, counts, -(1:3)) %>%
    left_join(ground_truth, by = c("subject", "locus")) %>%
    mutate(index = sub("est_counts.kallisto.", "", index),
           index = factor(index, levels = c("imgt", "pri")),
           prop_mapped = counts/true_counts) %>%
    rename(dist = dist.kallisto.imgt)

counts_star <- quant_data %>%
    select(subject, locus, dist.star.imgt, starts_with("est_counts.star")) %>%
    gather(index, counts, -(1:3)) %>%
    left_join(ground_truth, by = c("subject", "locus")) %>%
    mutate(index = sub("est_counts.star.", "", index),
           index = factor(index, levels = c("imgt", "pri")),
           prop_mapped = counts/true_counts) %>%
    rename(dist = dist.star.imgt)

sample_ids <- sprintf("sample_%02d", 1:50)

reads_lost_imgt_df <- 
    file.path("./PEreads_75bp/expression/star/mappings_2", 
              sample_ids, "reads_lost_to_other_genes_hla.tsv") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "subject") %>%
    complete(gene_from, gene_to, fill = list(perc = 0)) %>%
    filter(gene_from != gene_to) %>%
    group_by(gene_from, gene_to) %>%
    summarize(perc = mean(perc)) %>%
    ungroup() %>%
    filter(perc > 0) %>%
    mutate_at(vars(gene_from, gene_to), factor)

reads_gained_imgt_df <- 
    file.path("./PEreads_75bp/expression/star/mappings_2", 
              sample_ids, "reads_gained_from_other_genes_hla.tsv") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "subject") %>%
    complete(gene_to, gene_from, fill = list(perc = 0)) %>%
    filter(gene_from != gene_to) %>%
    group_by(gene_to, gene_from) %>%
    summarize(perc = mean(perc)) %>%
    ungroup() %>%
    filter(perc > 0) %>%
    mutate_at(vars(gene_to, gene_from), factor)

reads_lost_pri_df <- 
    file.path("./PEreads_75bp/expression/star/mappings_PRI", 
              sample_ids, "reads_lost_to_other_genes_hla.tsv") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "subject") %>%
    complete(gene_from, gene_to, fill = list(perc = 0)) %>%
    filter(gene_from != gene_to) %>%
    group_by(gene_from, gene_to) %>%
    summarize(perc = mean(perc)) %>%
    ungroup() %>%
    filter(perc > 0) %>%
    mutate_at(vars(gene_from, gene_to), factor)

reads_gained_pri_df <- 
    file.path("./PEreads_75bp/expression/star/mappings_PRI", 
              sample_ids, "reads_gained_from_other_genes_hla.tsv") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "subject") %>%
    complete(gene_to, gene_from, fill = list(perc = 0)) %>%
    filter(gene_from != gene_to) %>%
    group_by(gene_to, gene_from) %>%
    summarize(perc = mean(perc)) %>%
    ungroup() %>%
    filter(perc > 0) %>%
    mutate_at(vars(gene_to, gene_from), factor)

# Plots
png("./plots/kallisto_prop_mapped.png", width = 6, height = 4, units = "in", res = 200)
plot_dist(counts_kallisto)
dev.off()

png("./plots/star_prop_mapped.png", width = 6, height = 4, units = "in", res = 200)
plot_dist(counts_star)
dev.off()

png("./plots/kallisto_vs_star_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.kallisto.imgt", "tpm.star.imgt") +
    labs(x = "TPM (kallisto)", y = "TPM (STAR-Salmon)")
dev.off()

png("./plots/kallisto_vs_star_10pc.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "resid.kallisto.imgt", "resid.star.imgt") +
    labs(x = "PCA-corrected TPM (kallisto)", y = "PCA-corrected TPM (STAR-Salmon)")
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

png("./plots/kallisto_imgt_vs_PRI_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.kallisto.imgt", "tpm.kallisto.pri") +
    labs(x = "TPM (kallisto-IMGT)", y = "TPM (kallisto REF chromosomes)")
dev.off()

png("./plots/kallisto_imgt_vs_PRI_10pc.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "resid.kallisto.imgt", "resid.kallisto.pri") +
    labs(x = "PCA-corrected TPM (kallisto-IMGT)", 
         y = "PCA-corrected TPM (kallisto REF chromosomes)")
dev.off()

png("./plots/star_imgt_vs_PRI_TPM.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.star.imgt", "tpm.star.pri") +
    labs(x = "TPM (STAR-IMGT)", y = "TPM (STAR REF chromosomes)")
dev.off()

png("./plots/star_imgt_vs_PRI_10pc.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "resid.star.imgt", "resid.star.pri") +
    labs(x = "PCA-corrected TPM (STAR-IMGT)", 
         y = "PCA-corrected TPM (STAR REF chromosomes)")
dev.off()

png("./plots/reads_lost_imgt.png", width = 10, height = 6, units = "in", res = 200)
ggplot(reads_lost_imgt_df, aes(gene_from, gene_to)) +
    geom_point(aes(size = perc)) +
    theme_bw() +
    labs(x = "gene from", y = "gene to", size = "average percentage")
dev.off()

png("./plots/reads_gained_imgt.png", width = 10, height = 6, units = "in", res = 200)
ggplot(reads_gained_imgt_df, aes(gene_to, gene_from)) +
    geom_point(aes(size = perc)) +
    theme_bw() +
    labs(x = "gene to", y = "gene from", size = "average percentage")
dev.off()

png("./plots/reads_lost_pri.png", width = 10, height = 6, units = "in", res = 200)
ggplot(reads_lost_pri_df, aes(gene_from, gene_to)) +
    geom_point(aes(size = perc)) +
    theme_bw() +
    labs(x = "gene from", y = "gene to", size = "average percentage")
dev.off()

png("./plots/reads_gained_pri.png", width = 10, height = 6, units = "in", res = 200)
ggplot(reads_gained_pri_df, aes(gene_to, gene_from)) +
    geom_point(aes(size = perc)) +
    theme_bw() +
    labs(x = "gene to", y = "gene from", size = "average percentage")
dev.off()


devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

# Functions
 
plot_dist <- function(df) {

    ggplot(df, aes(dist, prop_mapped, color = index)) +
	    geom_point() +
	    geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
	              alpha = .8, size = 1.5) +
	    scale_x_continuous(labels = scales::percent,
	                       breaks = scales::pretty_breaks(n = 3)) +
    	scale_y_continuous(breaks = seq(0, 1.5, 0.5)) +
    	facet_wrap(~locus, scales = "free_x") +
    	scale_color_manual(
    	    values = c(hla = "#8491B4B2", 
    	               transcriptome = "#DC0000B2",
    	               genome_to_transc = "#7E6148FF", 
    	               genome = "#91D1C2B2"),
    	    labels = c(hla = "HLA-personalized",
    	               transcriptome = "Ref transcriptome",
    	               genome_to_transc = "Ref transcriptome v2",
    	               genome = "Ref genome")) +
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

hla_genes <- gencode_hla$gene_name

ground_truth <- read_tsv("./PEreads_75bp/data/phenotypes_counts_tpm.tsv") %>%
    filter(grepl("IMGT_(A|B|C|DPB1|DQA1|DQB1|DRB1)", Name)) %>%
    mutate(locus = sub("^IMGT_([^*]+).+$", "HLA-\\1", Name)) %>%
    group_by(subject, locus) %>%
    summarize(true_counts = sum(TrueCounts)) %>%
    ungroup()

kallisto_supplemented <- 
    "./PEreads_75bp/expression/kallisto/quantifications_2/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    mutate(allele = gsub("IMGT_", "", allele)) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), 
              tpm = sum(tpm)) %>%
    ungroup()

star_supplemented <- 
    "./PEreads_75bp/expression/star/supplemented/quantifications_2/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    mutate(allele = gsub("IMGT_", "", allele)) %>%
    left_join(allele_dist, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), 
              tpm = sum(tpm), 
              dist = mean(dist)) %>%
    ungroup()

star_transcriptome <- 
    "./PEreads_75bp/expression/star/transcriptome/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv()

star_genome <-
    "./PEreads_75bp/expression/star/genome/quantifications/compiled_hla_quants.tsv" %>%
    read_tsv() %>%
    rename(locus = gene)

star_genome_to_transcriptome <-
    "./PEreads_75bp/expression/star/genome_to_transcriptome/quantifications/processed_imgt_quants.tsv" %>% 
    read_tsv()

quant_data <-
    left_join(star_supplemented, kallisto_supplemented, 
              by = c("subject", "locus"), 
              suffix = c(".star.hla", ".kallisto.hla")) %>%
    left_join(star_transcriptome, by = c("subject", "locus")) %>%
    rename(est_counts.star.transcriptome = est_counts, 
           tpm.star.transcriptome = tpm) %>%
    left_join(star_genome_to_transcriptome, by = c("subject", "locus")) %>%
    rename(est_counts.star.genome_to_transc = est_counts,
           tpm.star.genome_to_transc = tpm) %>%
    left_join(star_genome, by = c("subject", "locus")) %>%
    rename(est_counts.star.genome = est_counts)

counts_star <- quant_data %>%
    select(subject, locus, dist, starts_with("est_counts.star")) %>%
    gather(index, counts, -(1:3)) %>%
    left_join(ground_truth, by = c("subject", "locus")) %>%
    mutate(index = sub("est_counts.star.", "", index),
           index = factor(index, levels = c("hla", "transcriptome", "genome_to_transc", "genome")),
           prop_mapped = counts/true_counts)

sample_ids <- sprintf("sample_%02d", 1:50)

alignments_to_diff_gene_supplemented <- 
    file.path("./PEreads_75bp/expression/star/supplemented/mappings_2", 
              sample_ids, "alignments_to_diff_gene_hla.tsv") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "subject") %>%
    complete(gene_from, gene_to, fill = list(perc = 0)) %>%
    filter(gene_from != gene_to) %>%
    group_by(gene_from, gene_to) %>%
    summarize(perc = mean(perc)) %>%
    ungroup() %>%
    filter(perc > 0)

alignments_to_diff_gene_transcriptome <- 
    file.path("./PEreads_75bp/expression/star/transcriptome/mappings", 
              sample_ids, "alignments_to_diff_gene_hla.tsv") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "subject") %>%
    complete(gene_from, gene_to, fill = list(perc = 0)) %>%
    filter(gene_from != gene_to) %>%
    group_by(gene_from, gene_to) %>%
    summarize(perc = mean(perc)) %>%
    ungroup() %>%
    filter(perc > 0)

alignments_to_diff_gene_df <- 
    list("HLA-personalized" = alignments_to_diff_gene_supplemented, 
         "Ref transcriptome" = alignments_to_diff_gene_transcriptome) %>%
    bind_rows(.id = "index") %>%
    mutate_at(vars(gene_from, gene_to), factor)

alignments_from_diff_gene_supplemented <- 
    file.path("./PEreads_75bp/expression/star/supplemented/mappings_2", 
              sample_ids, "alignments_from_diff_gene_hla.tsv") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "subject") %>%
    complete(gene_to, gene_from, fill = list(perc = 0)) %>%
    filter(gene_from != gene_to) %>%
    group_by(gene_to, gene_from) %>%
    summarize(perc = mean(perc)) %>%
    ungroup() %>%
    filter(perc > 0)

alignments_from_diff_gene_transcriptome <- 
    file.path("./PEreads_75bp/expression/star/transcriptome/mappings", 
              sample_ids, "alignments_to_diff_gene_hla.tsv") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "subject") %>%
    complete(gene_to, gene_from, fill = list(perc = 0)) %>%
    filter(gene_from != gene_to) %>%
    group_by(gene_to, gene_from) %>%
    summarize(perc = mean(perc)) %>%
    ungroup() %>%
    filter(perc > 0) 

alignments_from_diff_gene_df <- 
    list("HLA-personalized" = alignments_from_diff_gene_supplemented, 
         "Ref transcriptome" = alignments_from_diff_gene_transcriptome) %>%
    bind_rows(.id = "index") %>%
    filter(gene_to %in% gencode_hla$gene_name) %>%
    mutate_at(vars(gene_to, gene_from), factor)

# Plots
png("./plots/star_prop_mapped.png", width = 6, height = 4, units = "in", res = 200)
plot_dist(counts_star)
dev.off()

png("./plots/kallisto_vs_star.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.kallisto.hla", "tpm.star.hla") +
    labs(x = "TPM (kallisto)", y = "TPM (STAR-Salmon)")
dev.off()

png("./plots/star_HLA_vs_refTranscriptome.png", width = 10, height = 6, units = "in", res = 200)
scatter_plot(quant_data, "tpm.star.hla", "tpm.star.transcriptome") +
    labs(x = "TPM (HLA-personalized)", y = "TPM (Ref transcriptome)")
dev.off()

png("./plots/alignments_to_diff_gene.png", width = 12, height = 6, units = "in", res = 200)
ggplot(alignments_to_diff_gene_df, aes(gene_from, gene_to)) +
    geom_point(aes(size = perc)) +
    facet_wrap(~index) +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "gene from", y = "gene to", size = "average percentage")
dev.off()

png("./plots/alignments_from_diff_gene.png", width = 12, height = 6, units = "in", res = 200)
ggplot(alignments_from_diff_gene_df, aes(gene_to, gene_from)) +
    geom_point(aes(size = perc)) +
    facet_wrap(~index) +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "gene to", y = "gene from", size = "average percentage")
dev.off()

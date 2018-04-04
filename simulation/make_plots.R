devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)
 
# Proportion mapped vs divergence
 
pipelines <- c("Ref Transcriptome (ML)", "HLA-personalized (ML)",
               "Ref Transcriptome (Quasi)", "HLA-personalized (Quasi)", 
               "Conventional", "Ref Genome (Unique)") 

cols <- ggsci::pal_npg()(6) %>%
    setNames(pipelines)

allele_dist <- "~/hlaexpression/imgt_index_v2/distances_to_reference.tsv" %>%
    read_tsv()

hla_genes <- gencode_hla$gene_name

ground_truth <- read_tsv("./PEreads_75bp/data/phenotypes_trueCounts.tsv") %>%
    filter(grepl("IMGT_(A|B|C|DPB1|DQA1|DQB1|DRB1)", Name)) %>%
    mutate(locus = sub("^IMGT_([^*]+).+$", "HLA-\\1", Name)) %>%
    group_by(subject, locus) %>%
    summarize(true_counts = sum(TrueCounts)) %>%
    ungroup()

uniq_reads <-
    "./PEreads_75bp/expression/1-map_to_genome/quantifications_uniq/gw_quants.tsv" %>%
    read_tsv(col_names = FALSE) %>%
    inner_join(gencode_hla, c("X2" = "gene_id")) %>%
    select(subject = X1, locus = gene_name, est_counts = X3)

transcriptMap_ref <- 
    "./PEreads_75bp/expression/3-map_to_transcriptome/reference/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv()

transcriptMap_pers <- 
    "./PEreads_75bp/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
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
                             transcriptMap_pers = "HLA-personalized (ML)",
                             transcriptMap_ref = "Ref Transcriptome (ML)",
                             uniq = "Ref Genome (Unique)"),
           pipeline = factor(pipeline, levels = pipelines),
           prop_mapped = counts/true_counts)

png("./plots/prop_mapped_divergence.png", width = 6, height = 4, units = "in", res = 200)
ggplot(hla_counts_df, aes(dist, prop_mapped, color = pipeline)) +
    geom_point() +
    geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
              alpha = .8, size = 1.5) +
    scale_x_continuous(labels = scales::percent,
                       breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = seq(0, 1.5, 0.5)) +
    facet_wrap(~locus, scales = "free_x") +
    scale_color_manual(values = cols) +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          strip.text = element_text(size = 10, face = "bold"),
          legend.position = c(.75, .15)) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(x = "sequence divergence to the HLA reference allele", 
         y = "estimated counts / ground truth")
dev.off()

# Conventional vs Transcriptome vs Quasimapping

ground_truth_gw <- read_tsv("./PEreads_75bp/data/phenotypes_trueCounts.tsv") 

ground_truth_imgt <- ground_truth_gw %>%
    filter(grepl("IMGT_", Name)) %>%
    mutate(gene_name = sub("^IMGT_([^*]+).+$", "HLA-\\1", Name)) %>%
    left_join(select(gencode_pri_gene, gene_id, gene_name), by = "gene_name") %>%
    select(subject, gene_id, TrueCounts)

ground_truth_tx <- ground_truth_gw %>%
    filter(!grepl("IMGT", Name)) %>%
    left_join(gencode_pri_tx, by = c("Name" = "tx_id")) %>%
    select(subject, gene_id, TrueCounts)

ground_truth_gene <- bind_rows(ground_truth_tx, ground_truth_imgt) %>%
    group_by(subject, gene_id) %>%
    summarise(counts = sum(TrueCounts)) %>%
    ungroup()

conventional_all_genes <- 
    "./PEreads_75bp/expression/1-map_to_genome/quantifications.bed" %>%
    read_tsv() %>%
    select(gene_id = gid, sample_01:sample_50) %>%
    gather(sample, counts, -gene_id) %>%
    select(subject = sample, gene_id, counts)

refTranscript_all_genes <-
    "./PEreads_75bp/expression/3-map_to_transcriptome/reference/quantifications.bed" %>%
    read_tsv() %>%
    select(gene_id = gid, sample_01:sample_50) %>%
    gather(sample, counts, -gene_id) %>%
    select(subject = sample, gene_id, counts)

quasi_all_genes <-
    "./PEreads_75bp/expression/4-quasimapping/reference/quantifications.bed" %>%
    read_tsv() %>%
    select(gene_id = gid, sample_01:sample_50) %>%
    gather(sample, counts, -gene_id) %>%
    select(subject = sample, gene_id, counts)

all_genes_df <- list(conventional = conventional_all_genes, 
                     refTranscriptome = refTranscript_all_genes,
                     refTranscriptome_quasi = quasi_all_genes,
                     groundTruth = ground_truth_gene) %>%
    bind_rows(.id = "pipeline") %>%
    spread(pipeline, counts) %>%
    mutate_at(vars(conventional, refTranscriptome, refTranscriptome_quasi, groundTruth),
              ~ifelse(is.na(.), 0, .))

protein_coding_df <- all_genes_df %>%
    left_join(select(gencode_pri_gene, gene_id, gene_type), by = "gene_id") %>%
    filter(gene_type == "protein_coding") %>%
    select(-gene_type)

ref_mappings_df <-
    list("All genes" = all_genes_df, "Protein coding" = protein_coding_df) %>%
    bind_rows(.id = "gene_type")

cor_df <- ref_mappings_df %>%
    group_by(gene_type) %>%
    do(data.frame(r = cor(.$conventional, .$refTranscriptome),
                  p = cor(.$conventional, .$refTranscriptome, method = "spearman"))) %>%
    ungroup() %>%
    mutate_at(vars(r, p), ~round(., digits = 3)) %>%
    mutate(x = -2.5, yr = 4.5, yrho = 3.5)

png("./plots/ref_methods.png", width = 6, height = 3, units = "in", res = 200)
ggplot(ref_mappings_df, aes(log10(conventional), log10(refTranscriptome))) +
    geom_abline() +
    geom_point(size = 0.5) +
    facet_wrap(~gene_type) +
    geom_text(data = cor_df, aes(x = x, y = yr, 
                                 label = paste("r == ", r)), parse = TRUE) +
    geom_text(data = cor_df, aes(x = x, y = yrho, 
                                 label = paste("rho == ", p)), parse = TRUE) +
    theme_bw() +
    labs(x = expression(paste("-log"[10], "Counts (genome)")),
         y = expression(paste("-log"[10], "Counts (transcriptome)")))
dev.off()

ref_mappings_true <- ref_mappings_df %>%
    filter(gene_type == "All genes") %>%
    select(-gene_type) %>%
    gather(pipeline, counts, conventional, refTranscriptome, refTranscriptome_quasi)

cor_df_true <- ref_mappings_true %>%
    group_by(pipeline) %>%
    do(data.frame(r = cor(.$counts, .$groundTruth),
                  p = cor(.$counts, .$groundTruth, method = "spearman"))) %>%
    ungroup() %>%
    mutate_at(vars(r, p), ~round(., digits = 3)) %>%
    mutate(x = -2.5, yr = 4.5, yrho = 3.5)

png("./plots/ref_methods_groundTruth.png", width = 6, height = 3, units = "in", res = 200)
ggplot(ref_mappings_true, aes(log10(counts), log10(groundTruth))) +
    geom_abline() +
    geom_point(size = 0.5) +
    facet_wrap(~pipeline) +
    geom_text(data = cor_df_true, aes(x = x, y = yr, 
                                 label = paste("r == ", r)), parse = TRUE) +
    geom_text(data = cor_df_true, aes(x = x, y = yrho, 
                                 label = paste("rho == ", p)), parse = TRUE) +
    theme_bw() +
    labs(x = expression(paste("-log"[10], "Counts")),
         y = expression(paste("-log"[10], "Counts (ground truth)")))
dev.off()


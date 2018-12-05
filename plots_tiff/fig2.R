devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(scales)


allele_dist <- "../imgt_index/distance_to_ref/distances_to_reference.tsv" %>%
    read_tsv()

hla_genes <- gencode_hla$gene_name

hla_regex <- sub("HLA-", "", hla_genes) %>%
    paste(collapse = "|") %>%
    paste0("IMGT_(", ., ")")

ground_truth <- read_tsv("../simulation/data/phenotypes_trueCounts.tsv") %>%
    filter(grepl(hla_regex, Name)) %>%
    mutate(locus = sub("^IMGT_([^*]+).+$", "HLA-\\1", Name)) %>%
    group_by(subject, locus) %>%
    summarize(true_counts = sum(TrueCounts)) %>%
    ungroup()

uniq_reads <-
    "../simulation/expression/1-map_to_genome/quantifications_uniq/gw_quants.tsv" %>%
    read_tsv(col_names = FALSE) %>%
    inner_join(gencode_hla, c("X2" = "gene_id")) %>%
    select(subject = X1, locus = gene_name, est_counts = X3)

transcriptMap_ref <- 
    "../simulation/expression/3-map_to_transcriptome/reference/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv()

transcriptMap_pers <- 
    "../simulation/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% hla_genes) %>%
    mutate(allele = gsub("IMGT_", "", allele)) %>%
    left_join(allele_dist, by = c("locus", "allele")) %>%
    group_by(subject, locus) %>%
    summarize(est_counts = sum(est_counts), 
              tpm = sum(tpm), 
              dist = mean(dist)) %>%
    ungroup()

pipelines <- 
    c("HLApers (ML)", "Ref Transcriptome (ML)", "Ref Genome (Unique)") 

cols <- ggsci::pal_npg()(4) %>%
    .[c(1, 2, 4)] %>%
    setNames(pipelines)

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
    mutate(locus = sub("HLA-", "", locus),
           locus = factor(locus, levels = sub("HLA-", "", gencode_hla$gene_name)),
           pipeline = sub("^counts\\.", "", pipeline),
           pipeline = recode(pipeline, 
                             transcriptMap_pers = "HLApers (ML)",
                             transcriptMap_ref = "Ref Transcriptome (ML)",
                             uniq = "Ref Genome (Unique)"),
           pipeline = factor(pipeline, levels = pipelines),
           prop_mapped = counts/true_counts)


tiff("./plots/Fig2.tiff", width = 13.2, height = 8, units = "cm", res = 300)
ggplot(hla_counts_df, aes(dist, prop_mapped, color = pipeline)) +
    geom_line(stat = "smooth", method = "loess", span = 1, se = FALSE, 
              alpha = .4, size = 1) +
    geom_point(size = 1, alpha = .75) +
    scale_x_continuous(labels = function(x) round(x*100, 1),
                       breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5)) +
    facet_wrap(~locus, scales = "free_x") +
    scale_color_manual(values = cols) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme_bw() + 
    theme(text = element_text(family = "Times", size = 9), 
          strip.text = element_text(face = "bold")) +
    labs(x = "Sequence divergence to the HLA reference allele (%)", 
         y = expression(frac(Aligned~reads, Simulated~reads)))
dev.off()

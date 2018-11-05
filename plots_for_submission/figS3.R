library(tidyverse)
library(scales)
devtools::load_all("/home/vitor/Libraries/hlaseqlib/")

sample_info <- read_tsv("../geuvadis_reanalysis/replicates/data/sample_info.tsv")

replicates <- read_tsv("../geuvadis_reanalysis/replicates/3-map_to_transcriptome/hla_quantifications.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%    
    left_join(sample_info, by = c("subject" = "ena_id")) %>%
    mutate(allele = sub("IMGT_", "", allele),
           sample = "replicate") %>%
    select(subject = subject_id, sample, locus, allele, tpm)

final_set <- 
    "../geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%    
    left_join(sample_info, by = c("subject" = "ena_id")) %>%
    mutate(allele = sub("IMGT_", "", allele),
           sample = "finalset") %>%
    select(subject = subject_id, sample, locus, allele, tpm) %>%
    filter(subject %in% replicates$subject)

gene_df <- bind_rows(replicates, final_set) %>%
    group_by(subject, sample, locus) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup() %>%
    spread(sample, tpm)

cor_df <- gene_df %>%
    group_by(locus) %>%
    do(data.frame(r = cor(.$finalset, .$replicate),
                  p = cor(.$finalset, .$replicate, method = "spearman"),
                  x = min(.$finalset),
                  y = max(.$replicate))) %>%
    ungroup() %>%
    mutate_at(vars(r, p), ~round(., digits = 2)) %>%
    mutate(label = paste("r == ", r, "*','~rho ==", p))

tiff("./plots/S3_fig.tiff", width = 8, height = 6, units = "in", res = 300)
ggplot(gene_df, aes(finalset, replicate)) +
    geom_point() +
    scale_x_continuous(breaks = pretty_breaks(2)) +
    scale_y_continuous(breaks = pretty_breaks(2)) +
    facet_wrap(~locus, scales = "free") +
    geom_text(data = cor_df, aes(x, y, label = label), parse = TRUE, 
              hjust = "inward", vjust = "inward", size = 4) +
    theme_bw() 
dev.off()
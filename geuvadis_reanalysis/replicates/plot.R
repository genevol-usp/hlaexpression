library(tidyverse)
devtools::load_all("/home/vitor/Libraries/hlaseqlib/")

sample_info <- read_tsv("./data/sample_info.tsv")

replicates <- read_tsv("./3-map_to_transcriptome/hla_quantifications.tsv") %>%
    filter(locus %in% gencode_hla$gene_name) %>%    
    left_join(sample_info, by = c("subject" = "ena_id")) %>%
    mutate(allele = sub("IMGT_", "", allele),
           sample = "replicate") %>%
    select(subject = subject_id, sample, locus, allele, tpm)

final_set <- 
    "../expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%    
    left_join(sample_info, by = c("subject" = "ena_id")) %>%
    mutate(allele = sub("IMGT_", "", allele),
           sample = "finalset") %>%
    select(subject = subject_id, sample, locus, allele, tpm) %>%
    filter(subject %in% replicates$subject)
    

genos_rep <- select(replicates, subject, locus, allele)
genos_fin <- select(final_set, subject, locus, allele)

anti_join(genos_rep, genos_fin)

genos_rep %>% filter(subject == "NA11930", locus == "HLA-C")
genos_fin %>% filter(subject == "NA11930", locus == "HLA-C")
pag %>% filter(subject == "NA11930", locus == "C")

genos_rep %>% filter(subject == "HG00321", locus == "HLA-DQB1")
genos_fin %>% filter(subject == "HG00321", locus == "HLA-DQB1")
pag %>% filter(subject == "HG00321", locus == "DQB1")

genos_rep %>% filter(subject == "NA07357", locus == "HLA-C")
genos_fin %>% filter(subject == "NA07357", locus == "HLA-C")
pag %>% filter(subject == "NA07357", locus == "C")

genos_rep %>% filter(subject == "NA20773", locus == "HLA-DPB1")
genos_fin %>% filter(subject == "NA20773", locus == "HLA-DPB1")
final_set %>% filter(subject == "NA20773")
replicates %>% filter(subject == "NA20773")

genos_rep %>% filter(subject == "NA20538", locus == "HLA-C")
genos_fin %>% filter(subject == "NA20538", locus == "HLA-C")
pag %>% filter(subject == "NA20538", locus == "C")

genos_rep %>% filter(subject == "NA11830", locus == "HLA-B")
genos_fin %>% filter(subject == "NA11830", locus == "HLA-B")
pag %>% filter(subject == "NA11830", locus == "B")

genos_rep %>% filter(subject == "NA12842", locus == "HLA-A")
genos_fin %>% filter(subject == "NA12842", locus == "HLA-A")
pag %>% filter(subject == "NA12842", locus == "A")


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

png("./plots/replicates.png", width = 8, height = 6, units = "in", res = 300)
ggplot(gene_df, aes(finalset, replicate)) +
    geom_point() +
    scale_x_continuous(breaks = scales::pretty_breaks(2)) +
    scale_y_continuous(breaks = scales::pretty_breaks(2)) +
    facet_wrap(~locus, scales = "free") +
    geom_text(data = cor_df, aes(x, y, label = label), parse = TRUE, 
              hjust = "inward", vjust = "inward", size = 4) +
    theme_bw() 
dev.off()
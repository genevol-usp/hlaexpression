devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name

sample_ids <- sprintf("sample_%02d", 1:50)

alignments_to_diff_gene_supplemented <- 
    file.path("./PEreads_75bp/expression/star/old/supplemented/mappings_2", 
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
    file.path("./PEreads_75bp/expression/star/old/transcriptome/mappings", 
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
    file.path("./PEreads_75bp/expression/star/old/supplemented/mappings_2", 
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
    file.path("./PEreads_75bp/expression/star/old/transcriptome/mappings", 
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

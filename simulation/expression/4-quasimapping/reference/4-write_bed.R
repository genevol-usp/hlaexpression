devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

quant_dir <- "./quantifications"

gencode <- gencode_chr_tx %>%
    filter(chr %in% 1:22) %>%
    select(target_id = tx_id, gene_id, gene_name)

samples <- sprintf("sample_%02d", 1:50)

expression_df <- file.path(quant_dir, samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "subject") %>%
    select(subject, target_id = Name, counts = NumReads)

gene_df <- expression_df %>%
    inner_join(gencode, by = "target_id") %>%
    group_by(subject, gene_id) %>%
    summarise(counts = sum(counts)) %>%
    ungroup()

final_df <- gene_df %>%
    spread(subject, counts)

gene_bed <- inner_join(final_df, gencode_chr_gene, by = "gene_id") %>%
    filter(chr %in% 1:22) %>%
    mutate(chr = as.integer(chr), gid = gene_id) %>%
    select(`#chr` = chr, start, end, id = gene_id, gid, strd = strand, 
	   starts_with("sample")) %>%
    arrange(`#chr`, start)

write_tsv(gene_bed, "quantifications.bed")

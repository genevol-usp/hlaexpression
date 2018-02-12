devtools::load_all("~/hlaseqlib")
library(tidyverse)

quant_files <- sprintf("./quantifications/sample_%02d.gene.count.bed", 1:50)

out_df <- quant_files %>%
    map_df(. %>% read_tsv(col_types = "---c--i") %>%
    filter(gene %in% gencode_hla$gene_id) %>%
    gather(subject, est_counts, 2)) %>%
    left_join(select(gencode_hla, gene_id, gene_name), by = c("gene" = "gene_id")) %>%
    select(subject, gene = gene_name, est_counts) %>%
    arrange(subject, gene) 

write_tsv(out_df, "./quantifications/compiled_hla_quants.tsv")

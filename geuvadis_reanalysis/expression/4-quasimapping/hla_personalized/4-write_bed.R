devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

quant_dir <- "./quantifications"

gencode <- gencode_chr_tx %>%
    filter(chr %in% 1:22) %>%
    select(target_id = tx_id, gene_id, gene_name)

samples <- geuvadis_info %>%
    filter(kgp_phase3 == 1, pop != "YRI") %>%
    select(name, subject = ena_id)

expression_df <- file.path(quant_dir, samples$subject, "quant.sf") %>%
    setNames(samples$subject) %>%
    map_df(read_tsv, .id = "subject") %>%
    left_join(samples, by = "subject") %>%
    select(subject = name, target_id = Name, tpm = TPM)

autosomes_set <- inner_join(expression_df, gencode, by = c("target_id")) %>%
    select(subject, target_id, gene_id, tpm)

imgt_set <- expression_df %>%
    filter(grepl("^IMGT", target_id)) %>%
    mutate(gene_name = sub("^IMGT_([^\\*]+).+$", "HLA-\\1", target_id)) %>%
    inner_join(distinct(gencode, gene_id, gene_name), by = "gene_name") %>%
    select(subject, target_id, gene_id, tpm)

gene_set <- bind_rows(autosomes_set, imgt_set)

gene_df <- gene_set %>%
    group_by(subject, gene_id) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup() %>%
    complete(subject, gene_id, fill = list(tpm = 0))

expressed_genes <- gene_df %>%
    group_by(gene_id) %>%
    filter(mean(tpm>0.1) >= 0.5) %>%
    ungroup()
    
final_df <- expressed_genes %>%
    spread(subject, tpm)

gene_bed <- inner_join(final_df, gencode_chr_gene, by = "gene_id") %>%
    filter(chr %in% 1:22) %>%
    mutate(chr = as.integer(chr), gid = gene_id) %>%
    select(`#chr` = chr, start, end, id = gene_id, gid, strd = strand, 
	   starts_with("HG"), starts_with("NA")) %>%
    arrange(`#chr`, start)

write_tsv(gene_bed, "quantifications_expressed50%.bed")

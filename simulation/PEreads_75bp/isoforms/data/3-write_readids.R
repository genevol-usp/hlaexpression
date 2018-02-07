library(Biostrings)
library(tidyverse)
devtools::load_all("~/hlaseqlib")

hla_tx <- gencode_pri_tx %>%
    filter(gene_name %in% gencode_hla$gene_name) %>%
    select(gene_name, tx_id)

read_info <-
    system("zcat ./fastq/sample_01_1.fastq.gz | grep ^@read", intern = TRUE)

hla_read_info <- tibble(x = read_info) %>%
    mutate(tx_id = sub("^@read\\d+_([^ ]+) .+$", "\\1", x)) %>% 
    inner_join(hla_tx, by = "tx_id") %>%
    extract(x, c("readid", "m1", "m2"), "([^ ]+) mate1:([0-9-]+);mate2:([0-9-]+)") %>%
    separate(m1, c("m1.1", "m1.2"), sep = "-", convert = TRUE) %>%
    separate(m2, c("m2.1", "m2.2"), sep = "-", convert = TRUE) %>%
    rowwise() %>%
    mutate(pos = list(c(m1.1:m1.2, m2.1:m2.2))) %>%
    ungroup()

select(hla_read_info, gene_name, tx_id, readid, m1.1, m1.2, m2.1, m2.2) %>%
    mutate(readid = sub("^@", "", readid)) %>%
    arrange(gene_name, readid) %>%
    write_tsv("read_info.tsv")

hla_read_cov <- hla_read_info %>%
    select(gene_name, tx_id, pos) %>%
    unnest(pos) %>%
    count(gene_name, tx_id, pos)

write_tsv(hla_read_cov, "./hla_read_cov.tsv")

devtools::load_all("~/hlaseqlib")
library(tidyverse)

qtls <- read_qtltools("../3-conditional_analysis/conditional_70_all.txt.gz") %>%
    inner_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    filter(bwd_best == 1L) %>%
    select(locus = gene_name, rank, rsid = var_id, dist)

encode <- read_tsv("./hla.qtl.functional.elements.tsv")

regulome <- read_tsv("./regulomeDB_results.tsv")

out <- 
    left_join(encode, regulome, by = c("locus" = "gene", "rank", "rsid")) %>%
    left_join(qtls, by = c("locus", "rank", "rsid")) %>%
    select(locus, rank, rsid, dist, regulome_score = score, tfbs = tf, dhs, histone_marks) %>%
    mutate(dist = dist/1e3L)

write_tsv(out, "./results.tsv")

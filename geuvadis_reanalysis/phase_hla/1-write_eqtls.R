devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

eqtls <- 
    "../qtls/star/main_pipeline/supplemented/3-conditional_analysis/conditional_50_all.txt.gz" %>%
    read_qtltools() %>%
    filter(phen_id %in% gencode_hla$gene_id, bwd_best == 1L) %>%
    pull(var_id)

write_lines(eqtls, "./eqtl_rsID.txt")

devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)

hla_genes <- sort(gencode_hla$gene_name)

gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>% 
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_name, gene_id, start, end)

CRDactiv_bed <-
    read_delim("./LCL_CRDactivityQTL.txt", col_names = FALSE, delim = " ") %>%
    filter(X9 == "chr6", X20 < 0.05,
	   X10 >= min(gencode_hla_v19$start) - 1e6L,
	   X11 <= max(gencode_hla_v19$end) + 1e6L) %>%
    select(var_id = X8, phen_id = X1) %>%
    mutate(var_func = "CRDactiv") %>%
    unite(info, c("var_func", "phen_id"), sep = ":")

CRDstruct_bed <- 
    read_delim("./LCL_CRDstructureQTL.txt", col_names = FALSE, delim = " ") %>%
    filter(X9 == "chr6", X20 < 0.05,
	   X10 >= min(gencode_hla_v19$start) - 1e6L,
	   X11 <= max(gencode_hla_v19$end) + 1e6L) %>%
    select(var_id = X8, phen_id = X1) %>%
    mutate(var_func = "CRDstruc") %>%
    unite(info, c("var_func", "phen_id"), sep = ":")

rsmerge <- read_tsv("RsMergeArch.bcp.gz", col_names = FALSE) %>%
    select(rsHigh = X1, rsLow = X2, build_id = X3, rsCurrent = X7) %>%
    filter(build_id <= 149) %>%
    mutate_at(vars(rsHigh, rsLow, rsCurrent), function(x) paste0("rs", x)) %>%
    select(rsHigh, rsCurrent)

variants_df <- 
    bind_rows(list(CRDactiv_bed, CRDstruct_bed)) %>%
    left_join(rsmerge, by = c("var_id" = "rsHigh")) %>%
    mutate(var_id = ifelse(is.na(rsCurrent), var_id, rsCurrent)) %>%
    select(var_id, info)

write_tsv(variants_df, "./catalog.tsv", col_names = FALSE)

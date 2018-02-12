devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- sort(gencode_hla$gene_name)

gencode_hla_v19 <- "~/gencode_data/gencode.v19.annotation.gtf.gz" %>% 
    get_gencode_coords(feature = "gene") %>%
    filter(gene_name %in% hla_genes) %>%
    select(gene_name, gene_id, start, end)

#eQTL_bed <- 
#    read_delim("./LCL_eQTL.txt", col_names = FALSE, delim = " ") %>%
#    inner_join(select(gencode_hla_v19, gene_id, gene_name), by = c("X1" = "gene_id")) %>%
#    filter(X20 <= 0.05) %>%
#    select(var_id = X8, gene_name) %>%
#    mutate(var_func = "eQTL") %>%
#    unite(info, c("var_func", "gene_name"), sep = ":")

cQTL_bed <-
    read_delim("./LCL_cQTL.txt", col_names = FALSE, delim = " ") %>%
    filter(X2 == "chr6", 
	   X10 >= min(gencode_hla_v19$start) - 1e6,
	   X11 <= max(gencode_hla_v19$end) + 1e6,
	   grepl("^rs", X8)) %>%
    group_by(X1) %>%
    filter(X20 == min(X20)) %>%
    ungroup() %>%
    select(var_id = X8, phen_id = X1) %>%
    mutate(var_func = "cQTL") %>%
    unite(info, c("var_func", "phen_id"), sep = ":")

#CRDactiv_bed <-
#    read_delim("./LCL_CRDactivityQTL.txt", col_names = FALSE, delim = " ") %>%
#    filter(X9 == "chr6", X20 < 0.05,
#	   X10 >= min(gencode_hla_v19$start) - 1e6,
#	   X11 <= max(gencode_hla_v19$end) + 1e6) %>%
#    select(var_id = X8, phen_id = X1) %>%
#    mutate(var_func = "CRDactiv") %>%
#    unite(info, c("var_func", "phen_id"), sep = ":")
#
#CRDstruct_bed <- 
#    read_delim("./LCL_CRDstructureQTL.txt", col_names = FALSE, delim = " ") %>%
#    filter(X9 == "chr6", X20 < 0.05,
#	   X10 >= min(gencode_hla_v19$start) - 1e6,
#	   X11 <= max(gencode_hla_v19$end) + 1e6) %>%
#    select(var_id = X8, phen_id = X1) %>%
#    mutate(var_func = "CRDstruc") %>%
#    unite(info, c("var_func", "phen_id"), sep = ":")

#CRDexpression_bed <- 
#    read_delim("./LCL_CRDexpressionQTL.txt", col_names=FALSE, delim=" ") %>%
#    filter(X1 %in% gencode_hla_v19$gene_id, X14 == 1L) %>%
#    select(var_chr = X9, var_from = X10, var_to = X11, var_id = X8)

rsmerge <- read_tsv("RsMergeArch.bcp.gz", col_names = FALSE) %>%
    select(rsHigh = X1, rsLow = X2, build_id = X3, rsCurrent = X7) %>%
    mutate_at(vars(rsHigh, rsLow, rsCurrent), function(x) paste0("rs", x))

variants_df <- 
#    bind_rows(list(eQTL_bed, cQTL_bed, CRDactiv_bed, CRDstruct_bed)) %>%
    cQTL_bed %>%
    left_join(rsmerge, by = c("var_id" = "rsHigh")) %>%
    mutate(var_id = ifelse(is.na(rsCurrent), var_id, rsCurrent)) %>%
    select(var_id, info)

write_tsv(variants_df, "./variants.tsv", col_names = FALSE)

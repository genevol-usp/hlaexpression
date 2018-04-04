library(tidyverse)

dir.create("./covariates")

pcs_df <-  
    read_delim("./phenotypes_eur_pcs.pca", n_max = 100L, delim = " ") %>%
    mutate(SampleID = sub("^.+(PC\\d+)$", "\\1", SampleID)) %>%
    rename(id = SampleID)

out_basename <- "./covariates/covariates_pheno_"

for (pc in c(seq(5, 20, 5), seq(30, 100, 10))) {

    df_i <- filter(pcs_df, id %in% paste0("PC", seq_len(pc)))
    write_tsv(df_i, paste0(out_basename, pc, ".txt"))
}

devtools::load_all("~/hlaseqlib")
library(tidyverse)

make_genot_calls_df <- function(typings_df) {
    
    typings_df %>%
        mutate(subject = convert_ena_ids(as.character(subject)),
	       locus = sub("^HLA-", "", locus),
	       allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
	arrange(subject, locus, allele)
}

hla_genes <- gencode_hla$gene_name 

samples <- sprintf("sample_%02d", 1:50)

imgt_quants <- read_tsv("./quantifications_1/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = Name, 
	   est_counts = NumReads, tpm = TPM)

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

goldstd_genos <- read_tsv("../../../data/genos.tsv")

thresholds <- as.list(seq(0, .25, .05))
names(thresholds) <- seq(0, .25, .05)

typings <- 
  plyr::ldply(thresholds, 
	      function(th) hla_genotype_dt(imgt_quants, th),
	      .id = "th")

calls <- typings %>%
    filter(locus %in% hla_genes) %>%
    select(th, subject, locus, allele) %>%
    make_genot_calls_df() 

accuracies <- calls %>%
    split(.$th) %>%
    plyr::ldply(function(df) calc_genotyping_accuracy(df, goldstd_genos),
		.id = "th") %>%
    group_by(th) %>%
    mutate(th_average = mean(accuracy)) %>%
    ungroup()
  
write_tsv(accuracies, "./genotyping_accuracies_1.tsv")

best_th <- accuracies %>%
  slice(which.max(th_average)) %>%
  pull(th) %>%
  as.character()

out_df <- filter(typings, th == best_th) %>% select(-th)

write_tsv(out_df, "./quantifications_1/processed_imgt_quants.tsv")

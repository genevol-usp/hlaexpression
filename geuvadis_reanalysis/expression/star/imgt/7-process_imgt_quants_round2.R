devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- sort(gencode_hla$gene_name) 

make_genot_calls_df <- function(typings_df) {
    
    typings_df %>%
        mutate(subject = convert_ena_ids(as.character(subject)),
	       locus = sub("^HLA-", "", locus),
	       allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
	arrange(subject, locus, allele)
}

samples <- geuvadis_info %>% 
    filter(kgp_phase3 == 1L & pop != "YRI") %>%
    pull(ena_id)

imgt_quants <- read_tsv("./quantifications_1/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = Name,
	   est_counts = NumReads, tpm = TPM)

missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

goldstd_genos <- mutate(pag, allele = hla_trimnames(allele, 3))

out_df <- hla_genotype_dt(imgt_quants, th = 0) %>%
    hla_apply_zigosity_threshold(th = 0.2)

calls <- out_df %>%
    filter(locus %in% hla_genes) %>%
    select(subject, locus, allele) %>%
    make_genot_calls_df()

accuracies <- calc_genotyping_accuracy(calls, goldstd_genos)

write_tsv(accuracies, "./genotyping_accuracies_2.tsv")

write_tsv(out_df, "./quantifications_2/processed_imgt_quants.tsv")

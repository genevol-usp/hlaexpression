devtools::load_all("~/hlaseqlib")
library(tidyverse)

hla_genes <- gencode_hla$gene_name 

#samples <- sprintf("sample_%02d", 1:50)
samples <- "sample_01"

imgt_quants <- read_tsv("./quantifications_1/imgt_quants.tsv") %>%
    mutate(locus = imgt_to_gname(Name),
	   gene_id = gname_to_gid(locus)) %>%
    select(subject, locus, gene_id, allele = Name, 
	   est_counts = NumReads, tpm = TPM)

imgt_quants %>% filter(grepl("IMGT_A\\*01:01:01", allele))
    filter(locus %in% hla_genes) %>%
    group_by(locus) %>%
    top_n(20, est_counts) %>% 
    filter(locus == 'HLA-A') 
%>% 
    arrange(desc(est_counts)) %>%
    print(n = Inf)

index <- 
    "~/hlaexpression/imgt_index_v2/gencode.v25.PRI.IMGT.transcripts.fa" %>%
    Biostrings::readDNAStringSet()

index[c("ENST00000384703.1", "ENST00000384780.1")]

index_chr <-
    "~/gencode_data/gencode.v25.transcripts.fa" %>%
    Biostrings::readDNAStringSet()
names(index_chr) <- sub("^([^\\|]+).*$", "\\1", names(index_chr))


index_chr[c("ENST00000384703.1", "ENST00000384780.1")]



ambig <- read_tsv("./quantifications_1/sample_01/aux_info/ambig_info.tsv") %>%
    mutate(id = names(index))




missing_files <- samples[! samples %in% imgt_quants$subject]

if (length(missing_files) > 0L) {
    stop(paste("missing files:", paste(missing_files, collapse = " ")))
}

goldstd_genos <- read_tsv("../../data/genos.tsv")

thresholds <- as.list(seq(0, .25, .05))
names(thresholds) <- seq(0, .25, .05)

typings <- 
  plyr::ldply(thresholds, 
	      function(th) hla_genotype_dt(imgt_quants, th),
	      .id = "th")

calls <- typings %>%
    filter(locus %in% hla_genes) %>%
    select(th, subject, locus, allele) %>%
    mutate(allele = hla_trimnames(gsub("IMGT_", "", allele), 3)) %>%
    arrange(subject, locus, allele)

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

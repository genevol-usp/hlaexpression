devtools::load_all("~/Libraries/hlaseqlib")
library(tidyverse)

hla_genos <- 
    "../expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    select(subject, locus, allele) %>%
    mutate(subject = convert_ena_ids(subject),
	   allele = gsub("IMGT_", "", allele)) %>%
    arrange(subject, locus, allele)

hla_genos_recoded <- code_alleles(hla_genos)

write_tsv(hla_genos_recoded, "./codes-phase.inp")

snp_genos <- read_tsv("./eqtl_snps.vcf", comment = "##") %>%
    select(ID, POS, starts_with("HG"), starts_with("NA")) %>%
    gather(subject, geno, -(ID:POS)) %>%
    filter(subject %in% hla_genos$subject) %>%
    mutate(marker_type = "S") %>%
    select(subject, locus = ID, pos = POS, marker_type, allele = geno)

genos_df <- hla_genos_recoded %>%
    group_by(subject, locus) %>%
    summarise(allele = paste(code, collapse = "|")) %>%
    ungroup() %>%
    left_join(gencode_hla, by = c("locus" = "gene_name")) %>%
    mutate(marker_type = "M") %>%
    select(subject, locus, pos = start, marker_type, allele) %>%
    bind_rows(snp_genos) %>%
    arrange(subject, pos) %>%
    mutate(ix = "1|2") %>%
    separate_rows(ix, allele, sep = "\\|")

var_df <- genos_df %>%
    distinct(pos, locus, marker_type) %>%
    arrange(pos) %>%
    mutate(known = ifelse(marker_type == "M", "*", 0))

write_tsv(var_df, "./input_loci.tsv")

n_ind <- n_distinct(genos_df$subject)
n_loci <- nrow(var_df)
P <- c("P", var_df$pos) %>% paste(collapse = " ")
m_types <- var_df$marker_type %>% paste(collapse = "")

haps <- genos_df %>%
    select(subject, pos, ix, allele) %>%
    arrange(subject, pos, ix) %>%
    spread(pos, allele) %>%
    select(-ix)

file.create("phase.inp")
file.create("phase.known")
write(n_ind, "phase.inp", append = TRUE)
write(n_loci, "phase.inp", append = TRUE)
write(P, "phase.inp", append = TRUE)
write(m_types, "phase.inp", append = TRUE)

haps %>%
plyr::d_ply(~subject, 
	    function(x) {
	      write(paste0("#", unique(x$subject)), "phase.inp", append = TRUE)
	      write(unlist(x[1, -1]), "phase.inp", ncolumns = n_loci, append = TRUE)
	      write(unlist(x[2, -1]), "phase.inp", ncolumns = n_loci, append = TRUE)
	      write(paste(var_df$known, collapse = ""), "phase.known", append = TRUE) 
	    })

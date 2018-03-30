devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

dist_to_ref <- "../../imgt_index_v2/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)

hla_dist <- 
    "../expression/3-map_to_transcriptome/hla_personalized/quantifications/processed_imgt_quants.tsv" %>%
    read_tsv() %>%
    filter(locus %in% gencode_hla$gene_name) %>%
    left_join(select(geuvadis_info, name, ena_id), by = c("subject" = "ena_id")) %>%
    select(subject = name, gene_name = locus, allele) %>%
    mutate(allele = sub("IMGT_", "", allele)) %>%
    arrange(subject, gene_name, allele) %>%
    left_join(dist_to_ref, by = "allele") %>%
    group_by(subject, gene_name) %>%
    summarize(dist = mean(dist)) %>%
    ungroup()

phen_best_imgt <- 
    "../eqtl_mapping/transcriptomemapping/hla_personalized/1-phenotypes/phenotypes_eur_60.bed.gz" %>%
    read_tsv() %>% 
    inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
    select(gene_name, HG00096:NA20828) %>%
    gather(subject, resid, -gene_name) %>%
    select(subject, gene_name, resid)

phen_best_conventional <- 
    "../eqtl_mapping/genomemapping/1-phenotypes/phenotypes_eur_70.bed.gz" %>%
    read_tsv() %>% 
    inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
    select(gene_name, HG00096:NA20828) %>%
    gather(subject, resid, -gene_name) %>%
    select(subject, gene_name, resid)

expression_df <- 
    list(personalized = phen_best_imgt, conventional = phen_best_conventional) %>%
    bind_rows(.id = "index")

eqtl_df <- read_tsv("./best_eqtl.tsv") %>%
    left_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(index, gene_name, rsid = var_id, slope = bwd_slope)

eqtl_info <- read_tsv("./best_eqtl_snps.vcf", comment = "##") %>%
    select(-`#CHROM`, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    gather(subject, genotype, -ID) %>%
    inner_join(eqtl_df, by = c("ID" = "rsid")) %>%
    select(index, subject, gene_name, rsid = ID, genotype, slope) %>%
    separate(genotype, c("1", "2"), sep = "\\|") %>%
    gather(hap, allele, `1`:`2`) %>%
    group_by(index, subject, gene_name, rsid, slope) %>%
    summarize(genotype = paste(sort(allele), collapse = "/"),
	      dosage = sum(as.integer(allele))) %>%
    ungroup()

df <-
    left_join(eqtl_info, expression_df, by = c("index", "subject", "gene_name")) %>%
    left_join(hla_dist, by = c("subject", "gene_name")) %>%
    mutate(index = recode(index, "imgt" = "HLA_personalized", "ref" = "Reference"))

cor_list <- select(df, index, gene_name, resid, dist, dosage) %>%
    split(list(.$index, .$gene_name)) %>%
    map(~select(., -index, -gene_name))

cor_df <-
    map(cor_list, ~tibble(r = cor(.$resid, .$dist),
			  r2 = cor(.$dist, .$dosage),
			  partial_r = ppcor::pcor(.)$estimate["resid", "dist"])) %>%
    bind_rows(.id = "id") %>%
    separate(id, c("index", "gene"), sep = "\\.") %>%
    rename("expression~divergence" = r, 
	   "divergence~genotype" = r2,
	   "expression~divergence_genotype" = partial_r)

write_tsv(df, "./integrated_data.tsv")
write_tsv(cor_df, "./correlations.tsv")

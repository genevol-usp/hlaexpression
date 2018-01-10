devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

dist_to_ref <- 
    "../../simulation/PEreads_75bp/data/distances_to_reference.tsv" %>%
    read_tsv() %>%
    select(-locus)

hla_dist <- 
    "../expression/star/imgt/quantifications_2/processed_imgt_quants.tsv" %>%
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

phen_best_imgt <- "../qtls/star/imgt/1-phenotypes/phenotypes_eur_60.bed.gz" %>%
    read_tsv() %>% 
    inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
    select(gene_name, HG00096:NA20828) %>%
    gather(subject, resid, -gene_name) %>%
    select(subject, gene_name, resid)

phen_best_pri <- "../qtls/star/pri/1-phenotypes/phenotypes_eur_60.bed.gz" %>%
    read_tsv() %>% 
    inner_join(gencode_hla, by = c("gid" = "gene_id")) %>%
    select(gene_name, HG00096:NA20828) %>%
    gather(subject, resid, -gene_name) %>%
    select(subject, gene_name, resid)

expression_df <- list(imgt = phen_best_imgt, ref = phen_best_pri) %>%
    bind_rows(.id = "index")

eqtl_df <- read_tsv("./best_eqtl.tsv") %>%
    left_join(gencode_hla, by = c("phen_id" = "gene_id")) %>%
    select(index, gene_name, variant = var_id)

eqtl_info <- read_tsv("./best_eqtl_snps.vcf", comment = "##") %>%
    select(-`#CHROM`, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    gather(subject, genotype, -ID) %>%
    inner_join(eqtl_df, by = c("ID" = "variant")) %>%
    select(index, subject, gene_name, rsid = ID, genotype) %>%
    separate(genotype, c("1", "2"), sep = "\\|") %>%
    gather(hap, allele, `1`:`2`) %>%
    group_by(index, subject, gene_name, rsid) %>%
    summarize(genotype = paste(sort(allele), collapse = "/"),
	      dosage = sum(as.integer(allele))) %>%
    ungroup()

df <-
    left_join(eqtl_info, expression_df, by = c("index", "subject", "gene_name")) %>%
    left_join(hla_dist, by = c("subject", "gene_name")) %>%
    mutate(index = recode(index, "imgt" = "HLA_personalized", "ref" = "Reference"))


cor_list <- select(df, index, gene_name, dosage, resid, dist) %>%
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

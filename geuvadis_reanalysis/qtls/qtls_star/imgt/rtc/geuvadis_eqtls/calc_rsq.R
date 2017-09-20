devtools::load_all("~/hlaseqlib")
library(tidyverse)

samples <- geuvadis_info %>% 
  filter(pop != "YRI", kgp_phase3 == 1) %>%
  pull(name)

geuvadis_qtls <- 
  read_tsv("./catalog.tsv", col_names = FALSE) %>%
  mutate(V2 = sub("^[^:]+:([^:]+):.+$", "\\1", V2)) %>%
  select(gene = V2, variant = V1)

star_qtls <- 
  read_tsv("../../../../plots/eqtls.tsv") %>%
  filter(rank == 0) %>%
  select(gene = phen_id, variant = var_id)

writeLines(unique(c(star_qtls$variant, geuvadis_qtls$variant)), 
	   "./star_and_geuvadis.txt")

system("./subset_vcf.sh")

qtls_df <- 
  bind_rows(list(STAR = star_qtls, GEUVADIS = geuvadis_qtls), .id = "source")

vcf <- 
  read_tsv("./hla_eqtls.vcf", comment = "##") %>%
  select(ID, samples)

haps_df <- 
  gather(vcf, subject, genotype, -1) %>% 
  left_join(qtls_df, by = c("ID" = "variant")) %>%
  select(gene, source, variant = ID, subject, genotype)

genotypes_df <-
  haps_df %>%
  separate(genotype, c("1", "2"), sep = "\\|") %>%
  gather(hap, allele, `1`:`2`) %>%
  arrange(subject, gene, source, hap)

haps_f <-
  genotypes_df %>%
  select(-variant) %>%
  spread(source, allele) %>%
  count(gene, GEUVADIS, STAR) %>%
  group_by(gene) %>%
  mutate(f_hap = n/sum(n)) %>%
  ungroup() %>%
  select(gene, allele_geuvadis = GEUVADIS, allele_star = STAR, f_hap)

alleles_f <- 
  genotypes_df %>%
  count(gene, source, allele) %>%
  group_by(gene, source) %>%
  mutate(f = n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  spread(source, f) %>%
  rename(f_geuvadis = GEUVADIS, f_star = STAR)

ld_df <-
  left_join(haps_f, 
	    select(alleles_f, gene, allele, f_geuvadis), 
	    by = c("gene", "allele_geuvadis" = "allele")) %>%
  left_join(select(alleles_f, gene, allele, f_star), 
	    by = c("gene", "allele_star" = "allele")) %>%
  mutate(D = f_hap - f_geuvadis*f_star,
	 rsq = D^2 / ((f_geuvadis * (1 - f_geuvadis)) * (f_star * (1 - f_star))))

write_tsv(ld_df, "./ld_star_geuvadis_bestSNPs.tsv")


devtools::load_all("/home/vitor/hlaseqlib")
library(tidyverse)

gene_df <- 
    "~/hlaexpression/geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/gene_quantifications.tsv" %>%
    read_tsv()

expressed_genes <- gene_df %>%
    group_by(gene_id) %>%
    filter(mean(tpm>0) >= 0.75) %>%
    ungroup()
    
final_df <- expressed_genes %>%
    spread(subject, tpm)

gene_bed <- inner_join(final_df, gencode_chr_gene, by = "gene_id") %>%
    mutate(chr = as.integer(chr), gid = gene_id) %>%
    select(`#chr` = chr, start, end, id = gene_id, gid, strd = strand, 
	   starts_with("HG"), starts_with("NA")) %>%
    arrange(`#chr`, start)

write_tsv(gene_bed, "phenotypes.bed")
system("bgzip phenotypes.bed && tabix -p bed phenotypes.bed.gz")

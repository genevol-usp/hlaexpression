library(data.table)
library(Biostrings)
devtools::load_all("/home/vitor/hlaseqlib")
  
abundances <- 
    fread("../../../../geuvadis_reanalysis/expression/star/pri/quantifications/ERR188021/quant.sf")

setkey(abundances, Name)

index <- readDNAStringSet("~/gencode_data/gencode.v25.PRI.transcripts.fa")
index <- index[width(index) >= 75]

tx <- intersect(names(index), abundances$Name)

index <- index[tx]
writeXStringSet(index, "./polyester_index.fa")

abundances <- abundances[tx]
abundances[, TrueCounts := round(NumReads/sum(NumReads) * 3e7)]
abundances[, TrueTPM := counts_to_tpm(TrueCounts, EffectiveLength)]

phenotypes_tpm_counts <- abundances[, .(Name, TrueCounts, TrueTPM)]
fwrite(phenotypes_tpm_counts, "./phenotypes_counts_tpm.tsv", sep = "\t", quote = FALSE)

phenotypes <- phenotypes_tpm_counts[, .(ERR188021 = TrueCounts)]
fwrite(phenotypes, "./phenotypes.tsv", sep = "\t", quote = FALSE)

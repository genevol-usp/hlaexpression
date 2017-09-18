devtools::load_all("~/hlaseqlib")
library(data.table)

gencode12 <- 
  setDT(get_gencode_coords("~/gencode_data/gencode.v12.annotation.gtf.gz")
      )[gene_name %in% paste0("HLA-", c("A", "B", "C", "DQA1", "DQB1", "DRB1"))
      ][, .(gene_id, gene_name)]

geuvadis <- fread("zcat < ../../../../../data/previous_qtls/geuvadis_eur_eqtls.txt.gz")
setnames(geuvadis, c(1, 4), c("snp_id", "gene_id"))

hla <- geuvadis[gencode12, on = .(gene_id)][grepl("^rs|^snp", snp_id)]

bed <- data.table(chr = paste0("chr", hla$V5),
		  start = as.integer(hla$V7) ,
		  end = as.integer(hla$V7) + 1L,
		  info = paste(hla$snp_id, hla$gene_name, round(hla$V12, 2), sep = ":"))

fwrite(bed, "./geuvadis_hlaQTLs.bed", col.names = FALSE, quote = FALSE, sep = "\t")

system("./liftover_geuvadis.sh")

bedhg38 <- 
  fread("./geuvadis_hlaQTLs_hg38.bed"
      )[, .(chr = V1, pos = V2, info = V4)
      ][, chr := as.integer(sub("chr", "", chr))]

vcf <- fread("zcat < ../../../../../data/1000G_sites/ALL.chr6_GRCh38_sites.20170504.vcf.gz",
	     skip = "#CHROM", select = 1:3)
setnames(vcf, c("chr", "pos", "id"))

m <- vcf[bedhg38, on = .(chr, pos), nomatch = 0L]

catalog <- m[, .(id, info)]

fwrite(catalog, "./catalog.tsv", col.names = FALSE, quote = FALSE, sep = "\t")

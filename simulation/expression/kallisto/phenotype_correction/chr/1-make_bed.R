library(data.table)
devtools::load_all("/home/vitor/hlaseqlib")

setDT(gencode_chr_gene)

gene_exp <- fread("../../kallisto_CHR_expressed90%.csv")
  
gene_exp <- melt(gene_exp, id = 1, measure = 2:ncol(gene_exp), 
		 variable.name = "gene_id", value.name = "tpm")

gene_exp <- dcast(gene_exp, gene_id ~ subject, value.var = "tpm")

gene_exp <- 
  gene_exp[gencode_chr_gene, on = .(gene_id), nomatch = 0L
         ][chr %in% 1:22
         ][, `:=`(gid = gene_id, gene_name = NULL)]

setcolorder(gene_exp, c("chr", "start", "end", "gene_id", "gid", "strand", 
			grep("^sample_", names(gene_exp), value = TRUE)))

setnames(gene_exp, c("chr", "gene_id", "strand"), c("#chr", "id", "strd"))

bed_out <- "./phenotypes.bed"
bed_out_gz <- paste0(bed_out, ".gz")
fwrite(gene_exp, bed_out, sep = "\t", quote = FALSE)
system(paste("bgzip", bed_out, "&& tabix -p bed", bed_out_gz))

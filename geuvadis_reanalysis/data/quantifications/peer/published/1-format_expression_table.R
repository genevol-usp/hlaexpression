devtools::load_all("~/hlaseqlib")
library(data.table)

setDT(geuvadis_info)

samples <- geuvadis_info[kgp_phase3 == 1L & pop != "YRI", assay_name]

GD462 <- fread("zcat < ../../GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")

GD660 <- 
  fread("zcat < ../../GD660.GeneQuantRPKM.txt.gz"
      )[Gene_Symbol %in% GD462$Gene_Symbol
      ][, c("Gene_Symbol", samples), with = FALSE]
      
setnames(GD660, sub("^([^\\.]+).*$", "\\1", names(GD660)))

GD660_l <- melt(GD660, id = 1, measure = 2:ncol(GD660), variable.name = "subject")

expressed_genes <- GD660_l[, .(mean(value > 0)), by = .(Gene_Symbol)][V1 >= 0.9]

GD660_l <- GD660_l[Gene_Symbol %in% expressed_genes$Gene_Symbol]

GD660_w <- dcast(GD660_l, subject ~ Gene_Symbol, value.var = "value")

fwrite(GD660_w, "./geuvadis_fpkms.csv")

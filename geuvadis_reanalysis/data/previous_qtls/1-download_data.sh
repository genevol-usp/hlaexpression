#!/bin/bash

#selected publications
#wget -O battle_eqtls.xls http://genome.cshlp.org/content/suppl/2013/11/04/gr.155192.113.DC1/Supplemental_Data_1.xls
#wget -O barreiro_eqtls.xlsx http://www.cell.com/cms/attachment/2073022834/2068642455/mmc4.xlsx
#wget ftp://jungle.unige.ch/SGX/QTLs/LCL_eQTL.txt

#regulomeDB
#wget http://www.regulomedb.org/downloads/RegulomeDB.dbSNP141.txt.gz

zcat RegulomeDB.dbSNP141.txt.gz |\
    awk '$1 == "chr6" && $2 > 28900000 && X2 < 34100000' |\
    grep -E "1a$|1b$|1c$|1d$|1e$|1f$" - > RegulomeDB.dbSNP141.mhc.qtl.txt

#haploreg
#wget http://archive.broadinstitute.org/mammals/haploreg/data/eqtls_v4.1.tsv.gz

# rs ID history table
#wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/data/organism_data/RsMergeArch.bcp.gz


#!/bin/bash

vcf=../../qtls/genotypes/eur_maf05.vcf.gz

zcat $vcf | grep -v "^#" |\
    awk '$1 == 6 && $2 > 28900000 && $2 < 34100000 {print $3}' > mhc_rsids_vcf.txt

#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools
vcf=/home/vitor/hlaexpression/geuvadis_reanalysis/data/genotypes/eur_maf05.vcf.gz

zcat ../3-conditional_analysis/conditional_70_all.txt.gz |\
    awk '$19 == 1 && $20 == 1' > ./conditional_eqtl_best_significant.txt 

awk '{print $8}' conditional_eqtl_best_significant.txt | sort | uniq > eqtl_IDs.txt

$bcftools view --include '%ID=@eqtl_IDs.txt' -O v $vcf |\
    grep -v "^#" |\
    awk '{print $1, $2, $3, $4, $5}' > eqtl_info.txt

Rscript integrate_eqtl_info.R

rm conditional_eqtl_best_significant.txt eqtl_IDs.txt eqtl_info.txt

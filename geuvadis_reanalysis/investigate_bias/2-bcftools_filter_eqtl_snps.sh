#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools
vcf=../data/genotypes/eur_maf05.vcf.gz

$bcftools view --include '%ID=@best_eqtl_rsID.txt' -O v -o best_eqtl_snps.vcf $vcf 

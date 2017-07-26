#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools
vcf=../../../qtls/genotypes/eur_maf05.vcf.gz

$bcftools view --include '%ID=@best_eqtls_rsID.txt' -O v -o best_eqtl_snps.vcf $vcf 

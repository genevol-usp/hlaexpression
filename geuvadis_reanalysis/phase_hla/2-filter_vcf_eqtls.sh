#!/bin/bash

bcftools=/home/vitor/Libraries/bcftools/bcftools
vcf=../data/genotypes/eur_maf05.vcf.gz

$bcftools view --include '%ID=@eqtl_rsID.txt' -O v -o eqtl_snps.vcf $vcf 

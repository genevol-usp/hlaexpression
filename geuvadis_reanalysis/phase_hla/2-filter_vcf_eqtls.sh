#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools
vcf=../qtls/genotypes/eur_maf05.vcf.gz

$bcftools view --include '%ID=@eqtl_rsID.txt' -O v -o eqtl_snps.vcf $vcf 

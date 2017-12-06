#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools
vcf=../../../../../data/1000G_genotypes/ALL.chr6_GRCh38.genotypes.20170504.vcf.gz

$bcftools view --include '%ID=@previous_qtls.txt' -o ./hla_eqtls.vcf -O v $vcf


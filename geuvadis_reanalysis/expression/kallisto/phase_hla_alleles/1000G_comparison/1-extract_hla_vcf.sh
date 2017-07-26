#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools
vcf=/home/vitor/hlaexpression/geuvadis_reanalysis/data/1000G_genotypes/ALL.chr6_GRCh38.genotypes.20170504.vcf.gz
regions=./hla_regions.bed

$bcftools view -R $regions --genotype ^miss -v snps -O v -o tmp.vcf $vcf 

$bcftools view --exclude '%ID=@../../../../qtls/genotypes/duplicate_ids.txt' -o hla_snps.vcf -O v tmp.vcf

rm tmp.vcf

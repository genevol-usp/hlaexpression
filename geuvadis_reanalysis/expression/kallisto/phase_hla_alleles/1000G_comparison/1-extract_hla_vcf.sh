#!/bin/bash

bcftools=/home/vitor/bcftools/bcftools
vcf=/home/vitor/hlaexpression/geuvadis_reanalysis/data/1000G_genotypes/ALL.chr6.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz
regions=./hla_regions.bed

$bcftools view -R $regions --genotype ^miss -v snps -O v -o hla_snps.vcf $vcf 

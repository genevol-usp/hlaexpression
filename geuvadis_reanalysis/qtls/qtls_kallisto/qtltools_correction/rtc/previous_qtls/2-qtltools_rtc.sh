#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

VCF=../../../../genotypes/eur_maf05.vcf.gz
BED=../../phenotypes/phenotypes_eur_60.bed.gz
COV=../../../../pca_genotypes/covariates_genos.txt
HOTSPOTS=../hotspots_hg38_fixed_distances.bed
CATALOG=./previous_eQTLs.tsv
COND=../../conditional_analysis/conditional_60_all.txt.gz
LOG=./log.txt
OUT=./rtc_results.txt

$QTLtools rtc --vcf $VCF --bed $BED --region 6:28000000-34000000\
  --cov $COV --hotspot $HOTSPOTS --gwas-cis $CATALOG $COND --normal\
  --conditional --log $LOG --out $OUT

#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools

VCF=../../../genotypes/eur_maf05.vcf.gz
SAMPLES=../../../genotypes/samples.eur
BED=../../../phenotypes/phenotypes_eur_60.bed.gz
COV=../../../covariates_genos.txt
HOTSPOTS=../hotspots_hg38_fixed_distances.bed
CATALOG=./catalog.tsv
COND=../../../permutations/conditional_60_all.txt.gz
LOG=./log.txt
OUT=./rtc_GeuvadisBestQTLs.txt

$QTLtools rtc --vcf $VCF --bed $BED --include-samples $SAMPLES\
  --region 6:28000000-34000000\
  --cov $COV --hotspot $HOTSPOTS --gwas-cis $CATALOG $COND --normal\
  --conditional --log $LOG --out $OUT

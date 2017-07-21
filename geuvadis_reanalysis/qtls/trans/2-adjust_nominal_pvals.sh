#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

SAMPLES=../genotypes/samples.eur
VCF=../genotypes/eur_maf05.vcf.gz
COV=../pca_genotypes/covariates_genos.txt
BED=../phenotypes/phenotypes_eur_60.bed.gz
NULL=./trans_nulldist.best.txt.gz
OUT=./trans_adjust

$QTLtools trans --vcf $VCF --cov $COV --bed $BED --adjust $NULL --normal\
  --threshold 0.1 --out $OUT

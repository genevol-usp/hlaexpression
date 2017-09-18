#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

SAMPLES=../../../genotypes/samples.eur
VCF=../../../genotypes/eur_maf05.vcf.gz
COV=../../../pca_genotypes/covariates_genos.txt
BED=../phenotypes/phenotypes_eur_60.bed.gz
OUT=./trans_nulldist

$QTLtools trans --vcf $VCF --cov $COV --bed $BED --sample 1000 --normal --out $OUT

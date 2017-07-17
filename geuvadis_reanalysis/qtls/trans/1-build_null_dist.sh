#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools

SAMPLES=../genotypes/samples.eur
VCF=../genotypes/eur_maf05.vcf.gz
COV=../covariates_genos.txt
BED=../phenotypes/phenotypes_eur_60.bed.gz
OUT=./trans_nulldist

$QTLtools trans --vcf $VCF --cov $COV --bed $BED --sample 1000 --normal --out $OUT

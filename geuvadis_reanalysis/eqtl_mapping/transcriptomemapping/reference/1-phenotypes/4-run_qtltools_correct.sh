#!/bin/bash

qtltools=/home/vitor/qtltools/bin/QTLtools

PC=60
BED=./phenotypes.bed.gz
COV=./covariates/covariates_pheno_$PC.txt
OUT=./phenotypes_$PC.bed

$qtltools correct --bed $BED --cov $COV --normal --out $OUT
bgzip $OUT && tabix -p bed $OUT.gz

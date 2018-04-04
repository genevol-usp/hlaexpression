#!/bin/bash

qtltools=/home/vitor/qtltools/bin/QTLtools

BED=./phenotypes.bed.gz

# 0 PCs
$qtltools correct --bed $BED --normal --out ./phenotypes_0.bed
bgzip ./phenotypes_0.bed && tabix -p bed ./phenotypes_0.bed.gz

# 5-100 PCs
for pcs in $(seq 10 10 100)
do
    COV=./covariates/covariates_pheno_$pcs.txt
    OUT=./phenotypes_$pcs.bed

    $qtltools correct --bed $BED --cov $COV --normal --out $OUT
    bgzip $OUT && tabix -p bed $OUT.gz
done

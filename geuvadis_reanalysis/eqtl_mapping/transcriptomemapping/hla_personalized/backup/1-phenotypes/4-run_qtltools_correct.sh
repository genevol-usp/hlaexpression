#!/bin/bash

qtltools=/home/vitor/qtltools/bin/QTLtools

BED=./phenotypes_eur.bed.gz

# 0 PCs
$qtltools correct --bed $BED --normal --out ./phenotypes_eur_0.bed
bgzip ./phenotypes_eur_0.bed && tabix -p bed ./phenotypes_eur_0.bed.gz

# 5-100 PCs
for pcs in $(seq 5 5 20; seq 30 10 100)
do
    COV=./covariates/covariates_pheno_$pcs.txt
    OUT=./phenotypes_eur_$pcs.bed

    $qtltools correct --bed $BED --cov $COV --normal --out $OUT
    bgzip $OUT && tabix -p bed $OUT.gz
done

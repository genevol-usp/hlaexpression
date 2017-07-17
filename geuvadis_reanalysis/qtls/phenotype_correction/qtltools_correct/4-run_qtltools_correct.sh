#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools

BED=./phenotypes_eur.bed.gz

for pcs in $(seq 0 5 100)
do
  COV=./covariates/covariates_pheno_$pcs.txt
  OUT=./corrected/phenotypes_eur_$pcs.bed
  
  $QTLtools correct --bed $BED --cov $COV --normal --out $OUT
  
  bgzip $OUT && tabix -p bed $OUT.gz
done

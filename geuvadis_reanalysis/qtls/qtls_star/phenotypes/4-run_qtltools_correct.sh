#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

BED=./phenotypes_eur.bed.gz

for pcs in $(seq 0 5 30; seq 40 10 100)
do
  COV=./covariates/covariates_pheno_$pcs.txt
  OUT=./phenotypes_eur_$pcs.bed
  
  $QTLtools correct --bed $BED --cov $COV --normal --out $OUT
  
  bgzip $OUT && tabix -p bed $OUT.gz
done

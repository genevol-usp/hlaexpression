#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64
parallel=/home/vitor/parallel

N_PCS=60
JOBS=60
SAMPLES=../genotypes/samples.eur
VCF=../genotypes/eur_maf05.vcf.gz
COV=../pca_genotypes/covariates_genos.txt
BED=../phenotype_correction/qtltools_correct/corrected/phenotypes_eur_$N_PCS.bed.gz
THR=../permutations/results/permutations_$N_PCS.thresholds.txt
OUT=./conditional_$N_PCS
LOG=./log.txt
CMD_FILE=./cmd_cond.txt

for j in $(seq 1 $JOBS)
do
  echo $QTLtools cis --vcf $VCF --bed $BED --cov $COV --mapping $THR \
    --include-samples $SAMPLES --normal --chunk $j $JOBS \
    --out $OUT\_$j.txt --log $LOG
done>$CMD_FILE

$parallel --gnu -j $JOBS :::: $CMD_FILE

cat $OUT\_*.txt | gzip -c > $OUT\_all.txt.gz
rm $OUT\_*.txt $CMD_FILE

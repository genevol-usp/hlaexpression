#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools
parallel=/home/vitor/parallel

K=55
JOBS=60
SAMPLES=../genotypes/samples.all
VCF=../genotypes/all_chrs_nodups_nomiss_biallelic_maf05_either_eur_yri.vcf.gz
COV=../covariates_all_individuals_geno.txt
BED=phenotypes/phenotypes_all_individuals_$K.bed.gz
THR=permutations/permutations_$K\_all.thresholds.txt
OUT=permutations/conditional_$K
LOG=log/log_conditional.txt
CMD_FILE=cmd_cond.txt

echo > $CMD_FILE

for j in $(seq 1 $JOBS)
do
  echo $QTLtools cis --vcf $VCF --bed $BED --cov $COV --mapping $THR \
    --include-samples $SAMPLES --normal --chunk $j $JOBS \
    --out $OUT\_$j\.txt --log $LOG >> $CMD_FILE
done

$parallel --gnu -j $JOBS :::: $CMD_FILE

cat $OUT\_*.txt | gzip -c > $OUT\_all.txt.gz
rm $OUT\_*.txt $CMD_FILE

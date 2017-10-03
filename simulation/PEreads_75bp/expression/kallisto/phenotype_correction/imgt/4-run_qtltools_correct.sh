#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

BED=./phenotypes.bed.gz
COV=./covariates_pheno_10.txt
OUT=./phenotypes_10.bed

$QTLtools correct --bed $BED --cov $COV --normal --out $OUT

bgzip $OUT

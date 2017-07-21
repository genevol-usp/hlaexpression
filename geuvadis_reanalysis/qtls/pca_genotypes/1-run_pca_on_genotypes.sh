#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

VCF=./genotypes/eur_maf05.vcf.gz
OUT=./pca/eur
LOG=./pca/eur.log

$QTLtools pca --vcf $VCF --distance 6000 --center --scale --out $OUT --log $LOG

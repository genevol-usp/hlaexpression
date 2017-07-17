#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools

VCF=./genotypes/eur_maf05.vcf.gz
OUT=./pca/eur
LOG=./pca/eur.log

$QTLtools pca --vcf $VCF --distance 6000 --center --scale --out $OUT --log $LOG

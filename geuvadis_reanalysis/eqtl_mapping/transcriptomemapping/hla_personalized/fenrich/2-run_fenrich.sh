#!/bin/bash

QTLtools=/home/vitor/QTLtools/QTLtools_1.1_Ubuntu16.04_x86_64

BEDsignif=./genes_significant.bed
BEDquantif=./genes_quantified.bed
TFBS=./tfbs_encode_hg38.bed
OUT=./enrichment_qtl_in_tf.txt

$QTLtools fenrich --qtl $BEDsignif --tss $BEDquantif --bed $TFBS --out $OUT

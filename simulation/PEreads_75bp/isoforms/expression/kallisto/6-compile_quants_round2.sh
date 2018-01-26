#!/bin/bash

quantDir=./quantifications_2
OUTimgt=$quantDir/imgt_quants.tsv
sample=ERR188021
file=$quantDir/$sample/abundance.tsv

awk 'NR == 1 {print "subject\t" $0}' $file > $OUTimgt
awk -v SUBJECT="$sample" 'FNR > 1 && $1 ~ /IMGT/ {print SUBJECT "\t" $0}' $file >> $OUTimgt

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $OUTimgt > $OUTimgt.tmp &&\
    mv $OUTimgt.tmp $OUTimgt

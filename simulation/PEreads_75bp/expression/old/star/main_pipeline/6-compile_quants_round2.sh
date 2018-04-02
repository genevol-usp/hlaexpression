#!/bin/bash

quantDir=./quantifications_final
samples=(sample_{01..50})
OUTimgt=$quantDir/imgt_quants.tsv

awk 'FNR==1 {print "subject\t" $0}' $quantDir/${samples[0]}/quant.sf > $OUTimgt

for id in "${samples[@]}"
do
    file=$quantDir/$id/quant.sf

    awk -v SUBJECT="$id" 'FNR > 1 && $1 ~ /IMGT/ {print SUBJECT "\t" $0}' $file >> $OUTimgt
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $OUTimgt > $OUTimgt.tmp &&\
    mv $OUTimgt.tmp $OUTimgt

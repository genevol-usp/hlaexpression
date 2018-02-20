#!/bin/bash

quantDir=quantifications
samples=(sample_{01..50})
out=$quantDir/compiled_gw_quants.tsv

mkdir -p $quantDir
rm -f $out; touch $out

for id in "${samples[@]}"
do
    file=mappings/${id}_ReadsPerGene.out.tab 

    awk -v SUBJECT="$id" 'FNR>4 {print SUBJECT "\t" $0}' $file >> $out
done

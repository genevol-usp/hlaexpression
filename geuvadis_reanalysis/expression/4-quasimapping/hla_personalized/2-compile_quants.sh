#!/bin/bash

readarray -t samples < ../../../data/sample_info/samples_phase3_ena_eur.txt
quantIMGT=./quantifications
outIMGT=$quantIMGT/imgt_quants.tsv

awk 'FNR==1 {print "subject\t" $0}' $quantIMGT/${samples[0]}/quant.sf > $outIMGT

for id in "${samples[@]}"
do
    file=$quantIMGT/$id/quant.sf

    awk -v SUBJECT="$id" 'FNR > 1 && $1 ~ /IMGT/ {print SUBJECT "\t" $0}' $file >> $outIMGT
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $outIMGT > $outIMGT.tmp &&\
    mv $outIMGT.tmp $outIMGT

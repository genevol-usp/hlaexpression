#!/bin/bash

quantDir=$1

readarray -t samples < ../../data/sample_info/samples_phase3_ena_eur.txt

out=$quantDir/genomewide_quants.tsv

awk '{print "subject\t" $0}' <(head -n1 $quantDir/${samples[0]}/quant.sf) > $out

for id in "${samples[@]}" 
do
    file=$quantDir/$id/quant.sf
    awk -v SUBJECT="$id" 'NR != 1 {print SUBJECT "\t" $0}' $file >> $out
done


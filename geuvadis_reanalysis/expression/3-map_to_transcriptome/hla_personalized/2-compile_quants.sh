#!/bin/bash

quantDir=./quantifications

readarray -t samples < /home/vitor/hlaexpression/geuvadis_reanalysis/data/sample_info/samples_phase3_ena_eur.txt

OUTgw=$quantDir/all_transcripts_quants.tsv
OUTimgt=$quantDir/imgt_quants.tsv

awk 'NR == 1 {print "subject\t" $0}' $quantDir/${samples[0]}/quant.sf > $OUTgw

cp $OUTgw $OUTimgt

for id in "${samples[@]}" 
do
    file=$quantDir/$id/quant.sf

    awk -v SUBJECT="$id" 'NR != 1 {print SUBJECT "\t" $0}' $file >> $OUTgw

    awk -v SUBJECT="$id" 'NR != 1 && $1 ~ /IMGT/ {print SUBJECT "\t" $0}' $file >> $OUTimgt
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $OUTimgt > $OUTimgt.tmp &&\
    mv $OUTimgt.tmp $OUTimgt

gzip -f $OUTgw

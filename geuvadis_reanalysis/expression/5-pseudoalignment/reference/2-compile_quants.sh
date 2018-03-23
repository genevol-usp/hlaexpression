#!/bin/bash

quantDir=./quantifications
hlaIDs=/home/vitor/hlaexpression/index_transcriptome/hla_ids_pri.txt

readarray -t samples < ../../../data/sample_info/samples_phase3_ena_eur.txt

OUTgw=$quantDir/all_transcripts_quants.tsv
OUTimgt=$quantDir/imgt_quants.tsv

awk 'FNR == 1 {print "subject\t" $0}' $quantDir/${samples[0]}/abundance.tsv > $OUTgw

cp $OUTgw $OUTimgt

for id in "${samples[@]}" 
do
    file=$quantDir/$id/abundance.tsv
    
    awk -v SUBJECT="$id" 'FNR > 1 {print SUBJECT "\t" $0}' $file >> $OUTgw

    grep -F -f $hlaIDs $file |\
	awk -v SUBJECT="$id" 'FNR > 1 {print SUBJECT "\t" $0}' >> $OUTimgt
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $OUTimgt > $OUTimgt.tmp &&\
    mv $OUTimgt.tmp $OUTimgt

gzip -f $OUTgw

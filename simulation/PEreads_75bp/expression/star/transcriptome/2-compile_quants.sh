#!/bin/bash

quantDir=./quantifications
hlaIDs=/home/vitor/hlaexpression/imgt_index/hla_ids_pri.txt
samples=$(printf "sample_%02d" $PBS_ARRAYID)
OUTimgt=$quantDir/imgt_quants.tsv

awk 'NR == 1 {print "subject\t" $0}' $quantDir/${samples[0]}/quant.sf > $OUTimgt

for id in "${samples[@]}" 
do
    file=$quantDir/$id/quant.sf
    
    grep -F -f $hlaIDs $file |\
	awk -v SUBJECT="$id" 'FNR > 1 {print SUBJECT "\t" $0}' >> $OUTimgt
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $OUTimgt > $OUTimgt.tmp &&\
    mv $OUTimgt.tmp $OUTimgt

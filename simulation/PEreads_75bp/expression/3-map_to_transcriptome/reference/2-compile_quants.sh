#!/bin/bash

samples=(sample_{01..50})
hlaIDs=/home/vitor/hlaexpression/index_transcriptome/hla_ids_pri.txt
quant=./quantifications
out=$quant/imgt_quants.tsv

awk 'FNR==1 {print "subject\t" $0}' $quant/${samples[0]}/quant.sf > $out

for id in "${samples[@]}"
do
    file=$quant/$id/quant.sf

    grep -F -f $hlaIDs $file |\
	awk -v SUBJECT="$id" 'FNR>1 {print SUBJECT "\t" $0}' >> $out
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $out > $out.tmp && mv $out.tmp $out

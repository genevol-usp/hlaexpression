#!/bin/bash

samples=(sample_{01..50})
hlaIDs=/home/vitor/hlaexpression/index_transcriptome/hla_ids_pri.txt
quant=./quantifications
out=$quant/imgt_quants.tsv
quantUniq=./quantifications_uniq
outUniq=$quantUniq/gw_quants.tsv

awk 'FNR==1 {print "subject\t" $0}' $quant/${samples[0]}/quant.sf > $out

if [ -f "$outUniq" ]; then
    rm $outUniq
fi

touch $outUniq

for id in "${samples[@]}"
do
    file=$quant/$id/quant.sf

    grep -F -f $hlaIDs $file |\
	awk -v SUBJECT="$id" 'FNR>1 {print SUBJECT "\t" $0}' >> $out

    fileUniq=$quantUniq/${id}_ReadsPerGene.out.tab

    awk -v SUBJECT="$id" 'FNR>4 {print SUBJECT "\t" $1 "\t" $2}' $fileUniq >> $outUniq
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $out > $out.tmp && mv $out.tmp $out

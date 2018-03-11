#!/bin/bash

readarray -t samples < /home/vitor/hlaexpression/geuvadis_reanalysis/data/sample_info/samples_phase3_ena_eur.txt
hlaIDs=/home/vitor/hlaexpression/index_transcriptome/hla_ids_pri.txt
quantMHC=./quantifications_MHC
outMHC=$quantMHC/imgt_quants.tsv
quantTransc=./quantifications_transcriptome
outTransc=$quantTransc/imgt_quants.tsv

awk 'FNR==1 {print "subject\t" $0}' $quantMHC/${samples[0]}/quant.sf > $outMHC
awk 'FNR==1 {print "subject\t" $0}' $quantTransc/${samples[0]}/quant.sf > $outTransc

for id in "${samples[@]}"
do
    fileMHC=$quantMHC/$id/quant.sf
    fileTransc=$quantTransc/$id/quant.sf

    awk -v SUBJECT="$id" 'FNR>1 && $1 ~ /IMGT/ {print SUBJECT "\t" $0}' $fileMHC >> $outMHC
    
    grep -F -f $hlaIDs $fileTransc |\
	awk -v SUBJECT="$id" 'FNR>1 {print SUBJECT "\t" $0}' >> $outTransc
done

awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $outMHC > $outMHC.tmp &&\
    mv $outMHC.tmp $outMHC
awk '{print $1 "\t" $2 "\t" $5 "\t" $6}' $outTransc > $outTransc.tmp &&\
    mv $outTransc.tmp $outTransc

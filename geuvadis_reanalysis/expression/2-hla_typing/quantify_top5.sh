#!/bin/bash

salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon

sample=$1
CPUS=6
indexDIR=./sample_indices/$sample
sample_hla=./sample_indices/hla_$sample.fa

mkdir -p $indexDIR

$salmon index -t $sample_hla -i $indexDIR --type quasi -k 31

fq1=../1-map_to_genome/mappings/mhc_fqs/${sample}_1.fq
fq2=../1-map_to_genome/mappings/mhc_fqs/${sample}_2.fq
outQuant=./quantifications_top5
out=$outQuant/$sample

mkdir -p $out

$salmon quant -i $indexDIR -l IU -1 $fq1 -2 $fq2 -o $out -p $CPUS \
    --writeMappings > ${out}/mappings.sam

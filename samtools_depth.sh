#!/bin/bash

samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1
dir=./nci/expression/quantifications_2
bam=$dir/$sample\_Aligned.out.bam
bed=./nci/expression/imgt.bed
header=$dir/$sample\_header.sam
imgtbam=$dir/$sample\_imgt.bam
imgtmapped=$dir/$sample\_imgt_mapped.bam
imgtsorted=$dir/$sample\_imgt_sorted.bam
cov=$dir/$sample\_hla.cov

$samtools view -H $bam > $header

$samtools view $bam | awk -F $'\t' '$3 ~ /IMGT/' | cat $header - | $samtools view -Sb - > $imgtbam

$samtools view -b -f 0x2 -F 0x100 $imgtbam > $imgtmapped

$samtools sort $imgtmapped -o $imgtsorted

$samtools depth -a -m 100000 -b $bed $imgtsorted > $cov

rm $header $imgtbam $imgtmapped

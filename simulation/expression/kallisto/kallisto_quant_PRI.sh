#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto
samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1
fastqR1=../../data/fastq/${sample}_1.fastq.gz
fastqR2=../../data/fastq/${sample}_2.fastq.gz
index=../../../imgt_index/kallisto/gencode.v25.PRI.transcripts.idx
outdir=./quantifications_PRI
sampledir=$outdir/$sample
bam=$sampledir/alignments.bam

mkdir -p $sampledir

$kallisto quant -i $index -t 1 -o $sampledir --bias --pseudobam \
    $fastqR1 $fastqR2 | $samtools view -Sb - > $bam

awk 'NR==1 {print $1"\t"$4"\t"$5}' $sampledir/abundance.tsv >\
    $sampledir/abundance_imgt.tsv

grep -F -f ../../../imgt_index/hla_ids_pri.txt $sampledir/abundance.tsv |\
    awk -F $"\t" '{print $1"\t"$4"\t"$5}' >> $sampledir/abundance_imgt.tsv

header=$sampledir/header.sam
imgtbam=$sampledir/imgt.bam

$samtools view -H $bam > $header

$samtools view -f 0x2 $bam |\
    grep -F -f ../../data/ids_to_filter.txt - |\
    cat $header - |\
    $samtools view -Sb - > $imgtbam

Rscript ./parse_alignments.R $sample $imgtbam $sampledir

rm $bam

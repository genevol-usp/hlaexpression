#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1
mismatch=$2

indexDIR=./indexCHR
fq1=../../data/fastq/$sample\_1.fq.gz
fq2=../../data/fastq/$sample\_2.fq.gz
outMap=./mappings_CHR_$mismatch
outQuant=./quantifications_CHR_$mismatch

mkdir -p $outMap
mkdir -p $outQuant

outPrefix=$outMap/$sample\_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax $mismatch\
  --outFilterMultimapScoreRange 0\
  --outFilterMultimapNmax 50\
  --winAnchorMultimapNmax 100\
  --alignIntronMax 1\
  --alignEndsType Local\
  --outSAMunmapped Within KeepPairs\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

bam=$outPrefix\Aligned.out.bam
fasta=/home/vitor/gencode_data/gencode.v25.transcripts.fa
out=$outQuant/$sample

$salmon quant -t $fasta -l IU -a $bam -o $out -p 6

rm $outPrefix*

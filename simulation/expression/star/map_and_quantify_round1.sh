#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1

indexDIR=./index
fq1=../../data/fastq/$sample\_1.fq.gz
fq2=../../data/fastq/$sample\_2.fq.gz
outPrefix=./mappings_1/$sample\_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 1\
  --outFilterMultimapScoreRange 0\
  --outFilterMultimapNmax 4000\
  --winAnchorMultimapNmax 8000\
  --alignIntronMax 1\
  --alignEndsType EndToEnd\
  --outSAMunmapped None\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

bam=./mappings_1/$sample\_Aligned.out.bam
fasta=../../../geuvadis_reanalysis/expression/kallisto/index/gencode.v26.CHR.IMGT.transcripts.fa
out=./quantifications_1/$sample

$salmon quant -t $fasta -l IU -a $bam -o $out -p 6

rm ./mappings_1/$sample\_*

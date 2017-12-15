#!/bin/bash

STAR=/home/vitor/STAR

sample=$1

indexDIR=../../../../imgt_index/star/indexGENOME
fq1=../../data/fastq/${sample}_1.fastq.gz
fq2=../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings_GENOME
outPrefix=${outMap}/${sample}_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterType BySJout\
  --outFilterMultimapNmax 20\
  --winAnchorMultimapNmax 50\
  --alignIntronMax 0\
  --alignEndsType Local\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype None\
  --quantMode GeneCounts\
  --outFileNamePrefix $outPrefix


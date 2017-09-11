#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1

indexDIR=./indexPRI
fq1=../../data/fastq/$sample\_1.fastq.gz
fq2=../../data/fastq/$sample\_2.fastq.gz
outMap=./mappings_PRI
outQuant=./quantifications_PRI

mkdir -p $outMap
mkdir -p $outQuant

outPrefix=$outMap/$sample\_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterMultimapNmax 100\
  --winAnchorMultimapNmax 200\
  --alignIntronMax 0\
  --alignEndsType Local\
  --outSAMunmapped Within KeepPairs\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.out.bam
fasta=/home/vitor/gencode_data/gencode.v25.PRI.transcripts.fa
out=$outQuant/$sample

$salmon quant -t $fasta -l IU -a $bam -o $out -p 6

header=${outPrefix}header.sam
imgtbam=${outPrefix}imgt.bam

$samtools view -H $bam > $header

$samtools view -f 0x2 -F 0x100 $bam |\
  LC_ALL=C grep -F -f ./ids_to_filter.txt |\
  cat $header - |\
  $samtools view -Sb - > $imgtbam

rm ${outPrefix}Aligned* ${outPrefix}Log* ${outPrefix}SJ* ${outPrefix}header.sam

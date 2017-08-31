#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1
indexDIR=./sample_indices/$sample
fasta=../../../imgt_index/gencode.v25.PRI.transcripts.noIMGT.fa
sample_hla=./sample_indices/hla_$sample.fa
sample_fa=./sample_indices/index_$sample.fa

mkdir -p $indexDIR

cat $fasta $sample_hla > $sample_fa

$STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $indexDIR\
  --genomeFastaFiles $sample_fa --genomeChrBinNbits 11 --genomeSAindexNbases 13

fq1=../../data/fastq/$sample\_1.fastq.gz
fq2=../../data/fastq/$sample\_2.fastq.gz
outMap=./mappings_2
outQuant=./quantifications_2
outPrefix=$outMap/$sample\_

$STAR --runMode alignReads --runThreadN 16 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 0\
  --outFilterMultimapNmax 50\
  --winAnchorMultimapNmax 100\
  --alignIntronMax 0\
  --alignEndsType Local\
  --outSAMunmapped Within KeepPairs\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

rm -r $indexDIR

bam=${outPrefix}Aligned.out.bam
out=$outQuant/$sample

$salmon quant -t $sample_fa -l IU -a $bam -o $out -p 16

rm $sample_fa $sample_hla $bam ${outPrefix}Log.out ${outPrefix}Log.progress.out

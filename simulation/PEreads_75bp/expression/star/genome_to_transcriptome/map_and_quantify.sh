#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

indexDIR=/home/vitor/hlaexpression/index_genome/index
fq1=../../../data/fastq/${sample}_1.fastq.gz
fq2=../../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings
outPrefix=${outMap}/${sample}_

$STAR --runMode alignReads --runThreadN 8 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outSAMtype None\
  --quantMode TranscriptomeSAM\
  --quantTranscriptomeBan Singleend\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.toTranscriptome.out.bam
fasta=/home/vitor/gencode_data/gencode.v25.PRI.transcripts.fa
out=./quantifications/$sample

if [ -d "$out" ]; then
    rm -r $out
fi

$salmon quant -t $fasta -l IU -a $bam -o $out -p 8 --seqBias --gcBias

rm $outPrefix*

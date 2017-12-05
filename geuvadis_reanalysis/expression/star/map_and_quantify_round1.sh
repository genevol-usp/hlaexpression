#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1

indexDIR=/home/vitor/hlaexpression/imgt_index/star/index
fq1=../../data/fastq/${sample}_1.fastq.gz
fq2=../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings_1
outQuant=./quantifications_1
outPrefix=$outMap/${sample}_

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 1\
  --outFilterMultimapScoreRange 0\
  --outFilterMultimapNmax 3000\
  --winAnchorMultimapNmax 6000\
  --alignIntronMax 0\
  --alignEndsType EndToEnd\
  --outSAMunmapped None\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.out.bam
fasta=/home/vitor/hlaexpression/imgt_index/gencode.v25.PRI.IMGT.transcripts.fa
out=$outQuant/$sample

if [ -d "$out" ]; then
    rmdir $out
fi

$salmon quant -t $fasta -l IU -a $bam -o $out -p 6

awk 'NR==1 || $1 ~ /IMGT/ {print $1"\t"$4"\t"$5}' $out/quant.sf >\
    $out/quant_imgt.sf

rm $outPrefix*

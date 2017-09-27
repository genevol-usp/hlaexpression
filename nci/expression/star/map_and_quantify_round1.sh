#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1

indexDIR=/home/vitor/hlaexpression/imgt_index/star/index
fastqs1=$(ls -v ../../data/fastq/$sample*R1_001.fastq.gz)
fastqs2=$(ls -v ../../data/fastq/$sample*R2_001.fastq.gz)
outMap=./mappings_1
outQuant=./quantifications_1
outPrefix=$outMap/${sample}_

if [ "${#fastqs1[@]}" == 1 ]; then
    fq1=$fastqs1
    fq2=$fastqs2
elif [ "${#fastqs1[@]}" == 2 ]; then
    fq1=$(echo $fastqs1 | paste -d, - -)
    fq2=$(echo $fastqs2 | paste -d, - -)
else
    echo "wrong number of fastq files"
fi

$STAR --runMode alignReads --runThreadN 12 --genomeDir $indexDIR\
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

$salmon quant -t $fasta -l IU -a $bam -o $out -p 12

rm $outPrefix*

awk 'NR==1 || $1 ~ /IMGT/ {print $1"\t"$4"\t"$5}' $out/quant.sf >\
    $out/quant_imgt.sf

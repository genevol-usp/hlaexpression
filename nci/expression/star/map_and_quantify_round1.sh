#!/bin/bash

function join { local IFS="$1"; shift; echo "$*"; }

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

indexDIR=/home/vitor/hlaexpression/imgt_index/star/index
fastqs1=(`ls -v ../../data/fastq/$sample*R1_001.fastq.gz`)
fastqs2=(`ls -v ../../data/fastq/$sample*R2_001.fastq.gz`)
outMap=./mappings_1
outQuant=./quantifications_1
outPrefix=$outMap/${sample}_

if [ "${#fastqs1[@]}" == 1 ]; then
    fq1=$fastqs1
    fq2=$fastqs2
elif [ "${#fastqs1[@]}" == 2 ]; then
    fq1=$(join , ${fastqs1[@]})
    fq2=$(join , ${fastqs2[@]})
else
    echo "wrong number of fastq files"
fi

$STAR --runMode alignReads --runThreadN 8 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 1\
  --outFilterMultimapScoreRange 0\
  --outFilterMultimapNmax 1500\
  --winAnchorMultimapNmax 3000\
  --alignIntronMax 0\
  --alignEndsType EndToEnd\
  --outSAMunmapped None\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.out.bam
fasta=/home/vitor/hlaexpression/imgt_index/gencode.v25.PRI.IMGT.transcripts.fa
out=$outQuant/$sample

$salmon quant -t $fasta -l IU -a $bam -o $out -p 8

rm $outPrefix*

awk 'NR==1 || $1 ~ /IMGT/ {print $1"\t"$4"\t"$5}' $out/quant.sf >\
    $out/quant_imgt.sf

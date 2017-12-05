#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon

sample=$1

indexDIR=../../../imgt_index/star/indexPRI
fq1=../../data/fastq/${sample}_1.fastq.gz
fq2=../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings_PRI
outPrefix=$outMap/${sample}_
outQuant=./quantifications_PRI

$STAR --runMode alignReads --runThreadN 4 --genomeDir $indexDIR\
  --readFilesIn $fq1 $fq2 --readFilesCommand zcat\
  --outFilterMismatchNmax 999\
  --outFilterMismatchNoverReadLmax 0.04\
  --outFilterMultimapScoreRange 1\
  --outFilterMultimapNmax 150\
  --winAnchorMultimapNmax 300\
  --alignIntronMax 0\
  --alignEndsType Local\
  --outSAMunmapped Within KeepPairs\
  --outSAMprimaryFlag AllBestScore\
  --outSAMtype BAM Unsorted\
  --outFileNamePrefix $outPrefix

bam=${outPrefix}Aligned.out.bam
fasta=/home/vitor/gencode_data/gencode.v25.PRI.transcripts.fa
out=$outQuant/$sample

if [ -d "$out" ]; then
    rmdir $out
fi

$salmon quant -t $fasta -l IU -a $bam -o $out -p 4 --seqBias --gcBias

awk 'NR==1 {print $1"\t"$4"\t"$5}' $out/quant.sf > $out/quant_imgt.sf

grep -F -f ../../../imgt_index/hla_ids_pri.txt $out/quant.sf |\
    awk -F $"\t" '{print $1"\t"$4"\t"$5}' >> $out/quant_imgt.sf

rm ${outPrefix}Aligned* ${outPrefix}Log* ${outPrefix}SJ*

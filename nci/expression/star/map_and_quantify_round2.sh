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

$STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $indexDIR\
    --genomeFastaFiles $sample_fa\
    --genomeChrBinNbits 11 --genomeSAindexNbases 13\
    --outFileNamePrefix ${indexDIR}_

fq1=$(ls -v ../../data/fastq/$sample*R1_001.fastq.gz | paste -d, - -)
fq2=$(ls -v ../../data/fastq/$sample*R2_001.fastq.gz | paste -d, - -)
outMap=./mappings_2
outQuant=./quantifications_2
outPrefix=$outMap/${sample}_

mkdir -p $outMap
mkdir -p $outQuant

$STAR --runMode alignReads --runThreadN 8 --genomeDir $indexDIR\
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

rm -r $indexDIR ${indexDIR}_Log.out

bam=${outPrefix}Aligned.out.bam
out=$outQuant/$sample

$salmon quant -t $sample_fa -l IU -a $bam -o $out -p 8

awk 'NR==1 || $1 ~ /IMGT/ {print $1"\t"$4"\t"$5}' $out/quant.sf >\
    $out/quant_imgt.sf

rm $sample_fa $sample_hla $outPrefix* 

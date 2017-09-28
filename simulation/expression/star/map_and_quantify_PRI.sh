#!/bin/bash

STAR=/home/vitor/STAR
salmon=/home/vitor/Salmon-0.8.2_linux_x86_64/bin/salmon
samtools=/home/vitor/samtools-1.3.1/samtools

sample=$1

indexDIR=../../../imgt_index/star/indexPRI
fq1=../../data/fastq/${sample}_1.fastq.gz
fq2=../../data/fastq/${sample}_2.fastq.gz
outMap=./mappings_PRI
outPrefix=$outMap/${sample}_
outQuant=./quantifications_PRI

$STAR --runMode alignReads --runThreadN 6 --genomeDir $indexDIR\
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

$salmon quant -t $fasta -l IU -a $bam -o $out -p 6

awk 'NR==1 {print $1"\t"$4"\t"$5}' $out/quant.sf > $out/quant_imgt.sf

grep -F -f ../../../imgt_index/hla_ids_pri.txt $out/quant.sf |\
    awk -F $"\t" '{print $1"\t"$4"\t"$5}' >> $out/quant_imgt.sf

header=${outPrefix}header.sam
sampledir=$outMap/$sample
imgtbam=$sampledir/imgt.bam

$samtools view -H $bam > $header

$samtools view -f 0x2 -F 0x100 $bam |\
    grep -F -f ../../data/ids_to_filter.txt - |\
    cat $header - |\
    $samtools view -Sb - > $imgtbam

Rscript ./parse_alignments.R $sample $imgtbam $sampledir

rm ${outPrefix}Aligned* ${outPrefix}Log* ${outPrefix}SJ* ${outPrefix}header.sam

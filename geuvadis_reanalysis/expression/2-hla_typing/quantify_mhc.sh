#!/bin/bash

STAR=/home/vitor/STAR
samtools=/home/vitor/samtools-1.3.1/samtools
seqtk=/home/vitor/seqtk/seqtk
salmon=/home/vitor/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

CPUS=8
indexDIR=/home/vitor/hlaexpression/index_mhc/star/index
fq1=../1-map_to_genome/mappings/mhc_fqs/${sample}_1.fq
fq2=../1-map_to_genome/mappings/mhc_fqs/${sample}_2.fq
outMap=./mappings
outPrefix=${outMap}/${sample}_

$STAR --runMode alignReads --runThreadN $CPUS --genomeDir $indexDIR\
    --readFilesIn $fq1 $fq2\
    --outFilterMismatchNmax 1\
    --outFilterMultimapScoreRange 0\
    --outFilterMultimapNmax 3000\
    --winAnchorMultimapNmax 6000\
    --alignEndsType EndToEnd\
    --outSAMprimaryFlag AllBestScore\
    --outSAMtype BAM Unsorted\
    --outFileNamePrefix ${outPrefix}

fasta=/home/vitor/hlaexpression/index_mhc/gencode.v25.MHC.IMGT.transcripts.fa
bam=${outPrefix}Aligned.out.bam
out=./quantifications_MHC/$sample

if [ -d "$out" ]; then
    rm -r $out
fi

$salmon quant -t $fasta -l IU -a $bam -o $out -p $CPUS

rm ${outPrefix}*

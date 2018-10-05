#!/bin/bash

STAR=/home/vitor/Libraries/STAR
samtools=/home/vitor/Libraries/samtools-1.3.1/samtools
seqtk=/home/vitor/Libraries/seqtk/seqtk
salmon=/home/vitor/Libraries/Salmon-latest_linux_x86_64/bin/salmon

sample=$1

CPUS=8
indexDIR=/home/vitor/hlaexpression/index_mhc/star/index
outMap=./mappings
runid=$(pwgen 5 1)
outPrefix=/scratch/genevol/users/vitor/${runid}_${sample}_

fq1=../1-map_to_genome/mappings/mhc_fqs/${sample}_1.fq
fq2=../1-map_to_genome/mappings/mhc_fqs/${sample}_2.fq

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

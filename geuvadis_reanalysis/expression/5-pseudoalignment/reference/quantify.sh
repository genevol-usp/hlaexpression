#!/bin/bash

kallisto=/home/vitor/kallisto_linux-v0.43.1/kallisto

sample=$1
index=/home/vitor/hlaexpression/index_transcriptome/kallisto/gencode.v25.PRI.transcripts.idx
fq1=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_1.fastq.gz
fq2=/home/vitor/hlaexpression/geuvadis_reanalysis/data/fastq/${sample}_2.fastq.gz
out=quantifications/$sample

$kallisto quant -i $index -t 1 -o $out --bias $fq1 $fq2

#/bin/bash

bwakit=/home/vitor/bwa.kit

index=hs38DH.fa
fq1=../../data/fastq/sample_01_1.fastq.gz
fq2=../../data/fastq/sample_01_2.fastq.gz

out=sample_01

$bwakit/run-bwamem -o $out -H $index $fq1 $fq2 | sh

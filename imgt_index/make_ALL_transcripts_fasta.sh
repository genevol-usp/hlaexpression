#!/bin/bash

GENOMEGZ=/home/vitor/gencode_data/GRCh38.p7.genome.fa.gz
GENOME=/home/vitor/gencode_data/GRCh38.p7.genome.fa

GTFGZ=/home/vitor/gencode_data/gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz
GTF=/home/vitor/gencode_data/gencode.v25.chr_patch_hapl_scaff.annotation.gtf

OUT=./gencode.v25.ALL

zcat $GENOMEGZ > $GENOME
zcat $GTFGZ > $GTF

rsem-prepare-reference $GENOME $OUT --gtf $GTF

rm $OUT\.grp $OUT\.ti $OUT\.chrlist $OUT\.seq $OUT\.idx.fa $OUT\.n2g.idx.fa
rm $GENOME $GTF

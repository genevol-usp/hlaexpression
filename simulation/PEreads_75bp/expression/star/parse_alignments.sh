#!/bin/bash

sample=$1
outMap=./mappings_2
outPrefix=$outMap/${sample}_
bam=${outPrefix}Aligned.out.bam
sampledir=$outMap/$sample
imgtbam=$sampledir/imgt.bam

Rscript ./parse_alignments.R $sample $imgtbam $sampledir


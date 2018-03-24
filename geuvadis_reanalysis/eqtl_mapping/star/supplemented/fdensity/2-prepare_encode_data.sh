#!/bin/bash

liftDIR=/home/vitor/liftOver
liftover=$liftDIR/liftOver
chain=$liftDIR/hg19ToHg38.over.chain

# Files obtained from O. Delaneau:
TF=/home/vitor/ENCODE/TranscriptionFactors/step5_merged/ALL.bed.gz

OUT=TF.ENCODE.hg38.bed

$liftover $TF $chain $OUT ./unmapped_TF

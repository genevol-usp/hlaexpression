#!/bin/bash

FASTQ_DIR=../../data/fastq
SAMPLE=$1

zcat $FASTQ_DIR/$SAMPLE\_1.fastq.gz > $SAMPLE\_1.fastq &
zcat $FASTQ_DIR/$SAMPLE\_2.fastq.gz > $SAMPLE\_2.fastq
wait

python hlaTX.py --input $SAMPLE --outdir results --nthreads 4 --mismatch 0

for f in $(find . -name "$SAMPLE*"); do rm $f; done 

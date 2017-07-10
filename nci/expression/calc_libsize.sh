#!/bin/bash

fastqdir=/raid/genevol/rnaseq/fastq_NCI
out=library_size_NCI.txt

echo -n "" > $out

for i in $(find $fastqdir -name "*_R1_001.fastq.gz")
do
  filename="${i##*/}"
  sample="${filename:0:22}"

  echo $sample >> $out
  nlines=$(zcat $i | wc -l)
  nreads=$((nlines / 4))
  echo $nreads >> $out
done

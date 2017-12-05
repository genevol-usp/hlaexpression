#!/bin/bash

out=./library_size_NCI.txt

echo -n "" > $out

for i in $(find ./fastq -name "*_R1_001.fastq.gz")
do
  filename="${i##*/}"
  nlines=$(zcat $i | wc -l)
  nreads=$((nlines / 4))

  echo $filename$'\t'$nreads >> $out
done

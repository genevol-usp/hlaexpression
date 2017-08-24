#!/bin/bash

for sample in `cat ./samples_phase3_ena.txt`
do
  n=$(zcat ../fastq/$sample\_1.fastq.gz | grep -c ^@$sample)
  echo $sample$'\t'$n
done>./library_sizes.txt

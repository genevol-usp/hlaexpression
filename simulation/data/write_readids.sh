#!/bin/bash

mkdir -p read_ids

for i in sample_{01..50}
do
  zcat ./fastq/$i\_1.fq.gz | grep "^@" | grep "IMGT" | cut -d' ' -f1 > ./read_ids/$i.txt
done

#!/bin/bash

samples=(sample_{01..50})
annot=./quantifications/annotations.tsv
out=./quantifications/compiled_gw_quants.tsv

cut -f1-6 ./quantifications/${samples[0]}.gene.count.bed > $annot

for id in "${samples[@]}"
do
    cut -f7 ./quantifications/$id.gene.count.bed > ./quantifications/counts.$id
done

paste $annot ./quantifications/counts.sample_{01..50} > $out
rm ./quantifications/counts.sample_{01..50}

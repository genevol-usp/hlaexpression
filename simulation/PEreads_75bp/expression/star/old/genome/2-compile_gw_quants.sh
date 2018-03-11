#!/bin/bash

samples=(sample_{01..50})
annot=./quantifications/annotations.tsv
out=./quantifications/compiled_gw_quants.tsv
outTPM=./quantifications/compiled_gw_tpm.tsv
uniqHash=6OtGql2Sl8k

cut -f1-6 ./quantifications/${samples[0]}.$uniqHash.gene.count.bed > $annot
cut -f1-6 ./quantifications/${samples[0]}.$uniqHash.gene.tpm.bed > $annot

for id in "${samples[@]}"
do
    cut -f7 ./quantifications/$id.$uniqHash.gene.count.bed > ./quantifications/counts.$id
    cut -f7 ./quantifications/$id.$uniqHash.gene.tpm.bed > ./quantifications/tpm.$id
done

paste $annot ./quantifications/counts.sample_{01..50} > $out
paste $annot ./quantifications/tpm.sample_{01..50} > $outTPM

rm ./quantifications/counts.sample_{01..50}
rm ./quantifications/tpm.sample_{01..50}

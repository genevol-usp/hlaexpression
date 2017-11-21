#!/bin/bash

#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredInputsV3.tab.gz

zcat wgEncodeRegTfbsClusteredV3.bed.gz |\
    awk '{print $1"\t"$2"\t"$3"\t"$4}' > tfbs.bed

#zcat wgEncodeRegTfbsClusteredV3.bed.gz |\
#    awk '$1 == "chr6" {print $1"\t"$2"\t"$3"\t"$4}' > chr6_tfbs.bed

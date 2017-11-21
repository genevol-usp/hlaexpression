#!/bin/bash

liftdir=/home/vitor/liftOver
liftover=$liftdir/liftOver
chain=$liftdir/hg19ToHg38.over.chain

$liftover ./tfbs.bed $chain ./tfbs_hg38.bed ./unmapped_to_hg38

#$liftover ./chr6_tfbs.bed $chain ./chr6_tfbs_hg38.bed ./chr6_unmapped_to_hg38

#!/bin/bash

liftdir=/home/vitor/liftOver
liftover=$liftdir/liftOver
chain=$liftdir/hg19ToHg38.over.chain

$liftover ./chr6_tfbs.bed $chain ./chr6_tfbs_hg38.bed ./unmapped_to_hg38

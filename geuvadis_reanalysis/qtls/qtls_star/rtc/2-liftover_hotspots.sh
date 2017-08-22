#!/bin/bash

liftdir=/home/vitor/liftOver
liftover=$liftdir/liftOver
chain=$liftdir/hg19ToHg38.over.chain

$liftover hotspots_b37_hg19.bed $chain hotspots_hg38.bed unmapped_to_hg38

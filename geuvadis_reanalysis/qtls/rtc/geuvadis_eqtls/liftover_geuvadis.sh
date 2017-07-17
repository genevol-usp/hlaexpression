#!/bin/bash

liftdir=/home/vitor/liftOver
liftover=$liftdir/liftOver
chain=$liftdir/hg19ToHg38.over.chain

$liftover ./geuvadis_hlaQTLs.bed $chain ./geuvadis_hlaQTLs_hg38.bed unmapped_to_hg38

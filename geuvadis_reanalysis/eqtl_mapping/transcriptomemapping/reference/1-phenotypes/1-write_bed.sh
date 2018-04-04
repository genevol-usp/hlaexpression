#!/bin/bash

BED=../../../../expression/3-map_to_transcriptome/reference/quantifications_expressed50%.bed
OUT=./phenotypes_eur.bed

cp $BED $OUT
bgzip $OUT && tabix -p bed $OUT.gz

#!/bin/bash

BED=../../../../../expression/star/main_pipeline/quantifications_expressed50%_ref.bed
OUT=./phenotypes_eur.bed

cp $BED $OUT
bgzip $OUT && tabix -p bed $OUT.gz

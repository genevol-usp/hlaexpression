#!/bin/bash

BED=../../../../expression/star/supplemented/quantifications_expressed50%.bed
OUT=./phenotypes_eur.bed

cp $BED $OUT
bgzip $OUT && tabix -p bed $OUT.gz

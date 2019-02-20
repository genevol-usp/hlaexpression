#!/bin/bash

mkdir -p ./log

qtltools=/home/vitor/Libraries/qtltools/bin/QTLtools

PC=60
datadir=/home/vitor/hlaexpression/geuvadis_reanalysis/data
SAMPLES=$datadir/genotypes/samples.eur
VCF=$datadir/genotypes/eur_maf05.vcf.gz
COV=$datadir/pca_genotypes/covariates_genos.txt
HOTSPOTS=$datadir/recomb_hotspots/hotspots_hg38_fixed_distances.bed
CATALOG=$datadir/previous_qtls/top_qtl_catalog.tsv
BED=../../1-map_eqtls/th_50/1-phenotypes/phenotypes_$PC.bed.gz
COND=../../2-conditional_analysis/conditional_${PC}_all.txt.gz
OUT=./rtc_results.txt

$qtltools rtc --vcf $VCF --bed $BED --include-samples $SAMPLES\
    --region 6:29900000-33100000 --D-prime-threshold 0.1\
    --cov $COV --hotspot $HOTSPOTS --gwas-cis $CATALOG $COND --normal\
    --conditional --out $OUT

#!/bin/bash

qtltools=/home/vitor/qtltools/bin/QTLtools

$qtltools pca --bed ./phenotypes_eur.bed.gz --out ./phenotypes_eur_pcs --center --scale

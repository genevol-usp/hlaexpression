#!/bin/bash

qtltools=/home/vitor/qtltools/bin/QTLtools

$qtltools pca --bed ./phenotypes.bed.gz --out ./phenotypes_pcs --center --scale

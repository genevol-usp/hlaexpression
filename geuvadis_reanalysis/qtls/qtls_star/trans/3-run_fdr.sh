#!/bin/bash

runFDR=/home/vitor/QTLtools/script/runFDR_atrans.R
best=./trans_adjust.best.txt.gz
hits=./trans_adjust.hits.txt.gz
out=./fdr_pvals.txt

Rscript $runFDR $best $hits 0.05 $out

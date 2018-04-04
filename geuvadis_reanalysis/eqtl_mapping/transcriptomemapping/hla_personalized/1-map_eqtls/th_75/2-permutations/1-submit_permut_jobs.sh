#!/bin/bash

mkdir -p ./log
mkdir -p ./results

for PC in $(seq 0 10 100)
do
    for CHUNK in $(seq 1 64)
    do 
	qsub -v PC=$PC,CHUNK=$CHUNK run_qtltools_perm.pbs
    done
done

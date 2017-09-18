#!/bin/bash

for PC in $(seq 0 5 20; seq 30 10 100)
do
    for CHUNK in $(seq 1 60)
    do 
	qsub -v PC=$PC,CHUNK=$CHUNK run_qtltools_perm.pbs
    done
done

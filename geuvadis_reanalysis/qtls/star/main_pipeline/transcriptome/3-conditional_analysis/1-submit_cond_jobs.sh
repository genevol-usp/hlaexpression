#!/bin/bash

mkdir -p ./log

for CHUNK in $(seq 1 60)
do
    qsub -v PC=70,CHUNK=$CHUNK run_qtltools_conditional.pbs
done

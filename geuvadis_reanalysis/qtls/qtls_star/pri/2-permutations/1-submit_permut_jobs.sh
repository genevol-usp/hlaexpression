#!/bin/bash

for PC in $(seq 0 5 20; seq 30 10 100)
do
    for CHUNK in $(seq 1 60)
    do 
	qsub -l nodes=1:ppn=1,mem=1gb,walltime=24:00:00 -q short -j oe\
	    -o /home/vitor/hlaexpression/geuvadis_reanalysis/qtls/qtls_star/imgt/2-permutations/log/$PBS_JOBID.log\
	    -v PC=$PC,CHUNK=$CHUNK 1-run_qtltools_perm.sh
    done
done

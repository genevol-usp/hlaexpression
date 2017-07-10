#!/bin/bash

JOBS=48
parallel=/home/vitor/parallel

$parallel --gnu -j $JOBS ./kallisto_pseudobam.sh {} :::: ./nci_samples.txt

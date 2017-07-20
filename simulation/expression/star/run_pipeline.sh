#!/bin/bash

JOBS=10
parallel=/home/vitor/parallel
samples=$(echo sample_{01..50})

$parallel --gnu -j $JOBS ./2-map_and_quantify.sh {} ::: $samples

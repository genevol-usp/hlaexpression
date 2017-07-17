#!/bin/bash

/home/vitor/parallel --gnu --jobs 21 Rscript peer.R {} ::: $(seq 0 5 100) 

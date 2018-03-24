#!/bin/bash

/home/vitor/parallel --gnu --jobs 13 Rscript peer.R {} ::: $(seq 0 5 20; seq 30 10 100) 

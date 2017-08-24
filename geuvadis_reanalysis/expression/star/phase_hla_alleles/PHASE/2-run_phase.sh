#!/bin/bash

PHASE=/home/vitor/phase/src/phase.2.1.1.source/PHASE

$PHASE -d1 -x5 -F0.001 phase.inp phase.out 1000 1 100 &> log.txt

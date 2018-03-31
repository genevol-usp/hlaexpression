#!/bin/bash

awk 'FNR>1||NR==1' results/results{1..1842} > results.all

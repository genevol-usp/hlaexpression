#!/bin/bash

mv 4.1-caveman.pbs.o * ./log/
mv 4.2-caveman.pbs.o * ./log/

awk 'FNR>1||NR==1' results/results{1..1839} > results.all

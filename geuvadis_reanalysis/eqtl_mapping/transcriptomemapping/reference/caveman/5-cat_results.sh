#!/bin/bash

awk 'FNR>1||NR==1' results/results{1..1789} > results.all

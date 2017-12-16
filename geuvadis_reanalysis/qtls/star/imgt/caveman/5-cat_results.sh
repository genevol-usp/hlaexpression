#!/bin/bash

awk 'FNR>1||NR==1' results{1..920} > results.all

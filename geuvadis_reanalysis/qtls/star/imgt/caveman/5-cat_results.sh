#!/bin/bash

awk 'FNR>1||NR==1' results* > results.all

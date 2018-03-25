#!/bin/bash

zcat ../1-phenotypes/phenotypes_eur_60.bed.gz | cut -f-4,7- > phenotypes.bed

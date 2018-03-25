#!/bin/bash

zcat ../1-phenotypes/phenotypes_eur_50.bed.gz | cut -f-4,7- > phenotypes.bed

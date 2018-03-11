Report
================

Genotyping
==========

-   Threshold of least expressed/most expressed = 0.05

| locus    |  accuracy (%)|
|:---------|-------------:|
| HLA-A    |           100|
| HLA-B    |           100|
| HLA-C    |           100|
| HLA-DPB1 |            93|
| HLA-DQA1 |           100|
| HLA-DQB1 |           100|
| HLA-DRB1 |           100|

There are 7 mistypings for HLA-DPB1. In all of them, the true allele was DPB1\*04:01:01, but a DPB1\*126:01:01 was called instead. These 2 alleles are different only at position 619 of the CDS (end of exon 3). We need to investigate this issue.

Expression
==========

Pipelines:
----------

-   Ref Genome: primary assembly of reference genome GRCh38; mapped with STAR --quantMode GeneCounts
-   Ref transcriptome: reads mapped to Ref Genome with STAR, which creates alignments to transcripts given Gencode annotations; quantified with Salmon
-   HLA-personalized: Rescue of MHC-mapped and unmapped reads from the Ref Genome, alignment to MHC transcriptome supplemented with IMGT to infer HLA types, and final round of quantification of the whole-transcriptome supplemented with IMGT using Salmon

<img src="./plots/star_prop_mapped.png" width="1200" />

Comparisons between indices and aligners
========================================

kallisto vs STAR-Salmon; HLA-diversity index
--------------------------------------------

<img src="./plots/kallisto_vs_star.png" width="2000" />

STAR-Salmon; HLA-personalized vs Reference transcriptome
--------------------------------------------------------

<img src="./plots/star_HLA_vs_refTranscriptome.png" width="2000" />

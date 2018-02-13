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

-   Ref Genome: primary assembly of reference genome GRCh38; mapped with STAR, quantified with QTLtools quan
-   Ref transcriptome: transcripts from Ref Genome; mapped with STAR, quantified with Salmon
-   Ref transcriptome v2: reads mapped to Ref Genome with STAR, which creates alignments to transcripts given Gencode annotations; quantified with Salmon
-   HLA-personalized: Ref transcriptome where HLA transcripts were replace by IMGT CDS sequences; mapped with STAR, quantified with Salmon

<img src="./plots/star_prop_mapped.png" width="1200" />

**Note 1** To check: there is something strange about the genomic mappings above.

**Note 2** I also tried an alternative to replace the inefficient and time-consuming 1st round of the HLA-personalized pipeline. It consists in:

1.  Map to the genome, and filter the reads mapping to HLA-A, -C, -B, -DRB1, -DQA1, -DQB1 and -DPB1, or unmapped.
2.  Map these reads to the IMGT portion only.
3.  Quantify with Salmon.
4.  Type the HLA.

It works OK if I only filter for reads mapping to the genes described in step 1. However, if I filter reads mapping to all IMGT loci, step 2 fails for some reason. I need time to explore more this issue.

If this alternative works, it would save a lot of time, and allow people to use already existing BAM files.

The relevance of including all IMGT loci is that we could type and have better quantifications for them all. Our current pipeline does that.

Quality assessment
==================

Percentage of simulated reads not aligned, aligned exclusively to a different reference, or aligned to original reference but as a multimap
-------------------------------------------------------------------------------------------------------------------------------------------

| gene\_from | reads (%)                           |  HLA-personalized|  Ref transcriptome|
|:-----------|:------------------------------------|-----------------:|------------------:|
| HLA-A      | mapped to original gene as multimap |              1.21|               9.99|
| HLA-A      | uniquely mapping                    |             98.78|              80.58|
| HLA-A      | unmapped                            |              0.00|               4.67|
| HLA-A      | mapped but not to original gene     |              0.01|               4.76|
| HLA-B      | mapped to original gene as multimap |              0.36|               4.04|
| HLA-B      | uniquely mapping                    |             99.64|              74.22|
| HLA-B      | unmapped                            |              0.00|               9.18|
| HLA-B      | mapped but not to original gene     |              0.00|              12.56|
| HLA-C      | mapped to original gene as multimap |              0.01|               2.84|
| HLA-C      | uniquely mapping                    |             99.99|              79.00|
| HLA-C      | unmapped                            |              0.00|               4.15|
| HLA-C      | mapped but not to original gene     |              0.00|              14.02|
| HLA-DPB1   | mapped to original gene as multimap |              0.00|               0.00|
| HLA-DPB1   | uniquely mapping                    |            100.00|              97.81|
| HLA-DPB1   | unmapped                            |              0.00|               2.19|
| HLA-DPB1   | mapped but not to original gene     |              0.00|               0.01|
| HLA-DQA1   | mapped to original gene as multimap |              0.58|               9.75|
| HLA-DQA1   | uniquely mapping                    |             99.41|              59.79|
| HLA-DQA1   | unmapped                            |              0.00|              17.90|
| HLA-DQA1   | mapped but not to original gene     |              0.00|              12.56|
| HLA-DQB1   | mapped to original gene as multimap |              0.01|               2.34|
| HLA-DQB1   | uniquely mapping                    |             99.98|              70.40|
| HLA-DQB1   | unmapped                            |              0.00|              18.02|
| HLA-DQB1   | mapped but not to original gene     |              0.00|               9.24|
| HLA-DRB1   | mapped to original gene as multimap |              0.45|               4.05|
| HLA-DRB1   | uniquely mapping                    |             99.54|              61.84|
| HLA-DRB1   | unmapped                            |              0.00|              12.04|
| HLA-DRB1   | mapped but not to original gene     |              0.01|              22.07|

To further understand the read "loss" or "gain" we examined, for each HLA gene, the percentage of alignments stratified by mapped or source gene, respectively.

The size of each point in the plots below represents the percentage of alignments involving that gene pair.

### Percentage of aligments involving reads lost by an HLA gene in x to some other gene in y:

<img src="./plots/alignments_to_diff_gene.png" width="2400" />

### Percentage of alignments involving reads gained by an HLA gene in x from some other gene in y:

<img src="./plots/alignments_from_diff_gene.png" width="2400" />

Comparisons between indices and aligners
========================================

kallisto vs STAR-Salmon; HLA-diversity index
--------------------------------------------

<img src="./plots/kallisto_vs_star.png" width="2000" />

STAR-Salmon; HLA-personalized vs Reference transcriptome
--------------------------------------------------------

<img src="./plots/star_HLA_vs_refTranscriptome.png" width="2000" />

Isoforms
========

The simulation above is idealized in the sense that reads are generated from the coding sequences and mapped back to them. However, real transcription occurs via mRNA transcripts, containing UTRs and different combinations of exons.

In order to evaluate the performance of the HLA-supplemented index to map the isoforms, we simulated reads from the hypothetical homozygote individual for the reference HLA alleles at the 7 HLA loci. Ground truth isoform counts were based on the estimated counts for the individual NA20508 quantified with the Ref transcriptome index (protein coding transcripts only).

We see a reduction in expression in the TPM estimate for HLA-DQA1.

| gene     |  tpm (true)|  tpm (cds)|  estimated/true (cds)|
|:---------|-----------:|----------:|---------------------:|
| HLA-A    |    1283.142|   1300.923|                 1.014|
| HLA-B    |    2318.563|   2497.648|                 1.077|
| HLA-C    |     897.038|    948.303|                 1.057|
| HLA-DPB1 |     364.772|    344.513|                 0.944|
| HLA-DQA1 |     519.296|    405.828|                 0.781|
| HLA-DQB1 |     306.224|    308.304|                 1.007|
| HLA-DRB1 |    1061.513|   1123.560|                 1.058|

I tested versions of the index (1) including the UTRs or (2) trimming the index to 750bp, but they do not solve the problem (data not shown).

The problem is likely caused by shorter isoforms. When they are abundant, there is a underestimation of the gene TPMs, because normalization is done by a longer length (CDS).

If I don't filter for protein coding transcripts only, I also see a reduction for DQB1.

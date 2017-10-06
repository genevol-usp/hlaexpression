Report
================

Typing accuracies
=================

\*Concordance: the proportion of the called alleles that are concordant with the Gourraud et al (2014) typings

| locus |   star|  kallisto|
|:------|------:|---------:|
| A     |  97.55|     97.06|
| B     |  98.04|     97.55|
| C     |  96.90|     95.92|
| DQB1  |  98.37|     98.53|
| DRB1  |  99.51|     99.51|

Expression estimates
====================

![](./expression/plots/expression_boxplot.png)

kallisto vs STAR-Salmon
-----------------------

### TPM

![](./expression/plots/star_vs_kallisto_TPM.png)

### PCA-corrected

![](./expression/plots/star_vs_kallisto_PCA.png)

HLA diversity vs reference chromosomes only
-------------------------------------------

### TPM

#### STAR

![](./expression/plots/star_imgt_vs_pri_TPM.png)

#### kallisto

![](./expression/plots/kallisto_imgt_vs_pri_TPM.png)

### PCA-corrected

#### STAR

![](./expression/plots/star_imgt_vs_pri_PCA.png)

#### kallisto

![](./expression/plots/kallisto_imgt_vs_pri_PCA.png)

Distribution of TPM values
--------------------------

![](./expression/plots/tpm_distributions.png)

ASE
---

### ASE by number of genotyping errors

\*Each point represents a heterozygous genotype in the intersect with Gourraud data.

![](./expression/plots/ase.png)

### ASE distribution

![](./expression/plots/ase_histogram.png)

Correlation of expression
-------------------------

![](./expression/plots/correlation_decrease.png)

### Among the HLA genes

![](./expression/plots/hlacorrelations.png)

### Between Class II genes and CIITA

![](./expression/plots/trans_activ_corrs.png)

### Between pairs of HLA genes on the same vs on different haplotypes

#### HLA-A vs HLA-B

![](./expression/plots/a_vs_b.png)

#### HLA-A vs HLA-C

![](./expression/plots/a_vs_c.png)

#### HLA-B vs HLA-C

![](./expression/plots/b_vs_c.png)

#### HLA-DQA1 vs HLA-DQB1

![](./expression/plots/dqa_vs_dqb.png)

#### HLA-DQA1 vs HLA-DRB1

![](./expression/plots/dqa_vs_drb.png)

eQTLs
=====

**All analyses were carried out using European individuals only**

PCA of genotypes
----------------

![](./qtls/plots/genotype_pca.png)

Number of eGenes according to index
-----------------------------------

![](./qtls/plots/n_of_egenes.png)

Distribution of eQTLs around the TSS
------------------------------------

### IMGT index

![](./qtls/plots/qtls_landscape_imgt.png)

### Reference transcriptome

![](./qtls/plots/qtls_landscape_pri.png)

RegulomeDB score
----------------

| gene\_name |  rank| index | var\_id     | score |
|:-----------|-----:|:------|:------------|:------|
| HLA-A      |     0| imgt  | rs16896724  | 5     |
| HLA-A      |     0| pri   | rs28724990  | 6     |
| HLA-A      |     1| imgt  | rs2517827   | 7     |
| HLA-A      |     1| pri   | rs2916818   | 5     |
| HLA-A      |     2| imgt  | rs7383537   | 6     |
| HLA-A      |     2| pri   | rs2248162   | 4     |
| HLA-B      |     0| imgt  | rs1265094   | 4     |
| HLA-B      |     0| pri   | rs1265159   | 4     |
| HLA-B      |     1| imgt  | rs9264803   | 6     |
| HLA-B      |     1| pri   | rs2844623   | 6     |
| HLA-B      |     2| pri   | rs72443738  | 5     |
| HLA-B      |     3| pri   | rs28380929  | 4     |
| HLA-C      |     0| imgt  | rs41561715  | 4     |
| HLA-C      |     0| pri   | rs9468965   | 6     |
| HLA-C      |     1| imgt  | rs2394941   | 5     |
| HLA-C      |     1| pri   | rs2074491   | 4     |
| HLA-C      |     2| imgt  | rs2249741   | 1f    |
| HLA-C      |     2| pri   | rs2394941   | 5     |
| HLA-DPB1   |     0| imgt  | rs9277538   | 4     |
| HLA-DPB1   |     0| pri   | rs9277538   | 4     |
| HLA-DPB1   |     1| imgt  | rs9296068   | 6     |
| HLA-DPB1   |     1| pri   | rs9296068   | 6     |
| HLA-DPB1   |     2| imgt  | rs688209    | 4     |
| HLA-DQA1   |     0| imgt  | rs9272316   | 6     |
| HLA-DQA1   |     0| pri   | rs144150162 | NA    |
| HLA-DQA1   |     1| imgt  | rs9271375   | 7     |
| HLA-DQA1   |     1| pri   | rs35407265  | 3a    |
| HLA-DQB1   |     0| imgt  | rs1770      | 1f    |
| HLA-DQB1   |     0| pri   | rs28746846  | 3a    |
| HLA-DQB1   |     1| imgt  | rs9272209   | 3a    |
| HLA-DQB1   |     1| pri   | rs9273595   | 5     |
| HLA-DQB1   |     2| imgt  | rs9269386   | 7     |
| HLA-DQB1   |     2| pri   | rs114969562 | 5     |
| HLA-DQB1   |     3| pri   | rs3134988   | 7     |
| HLA-DRB1   |     0| imgt  | rs3104412   | 6     |
| HLA-DRB1   |     0| pri   | rs9269749   | 5     |
| HLA-DRB1   |     1| imgt  | rs9270700   | 5     |
| HLA-DRB1   |     1| pri   | rs9271365   | 7     |
| HLA-DRB1   |     2| pri   | rs182670197 | 7     |

Entire region
-------------

### Class I

![](./qtls/plots/landscape_class1.png)

### Class II

![](./qtls/plots/landscape_class2.png)

HLA lineages
------------

![](./qtls/plots/lineage_and_effects.png)

### F-test: do the grouping into lineages explain the variance?

| locus | term    |   df|     sumsq|     meansq|           F|  p.value|
|:------|:--------|----:|---------:|----------:|-----------:|--------:|
| A     | lineage |   14|  11855313|   846808.0|   10.523137|        0|
| B     | lineage |   25|  13384476|   535379.1|    4.858555|        0|
| C     | lineage |   13|  12001846|   923218.9|   27.398400|        0|
| DPB1  | lineage |   21|   2505297|   119299.8|    7.949565|        0|
| DQA1  | lineage |    5|  18367822|  3673564.4|  107.820537|        0|
| DQB1  | lineage |    4|   1491870|   372967.6|   37.588688|        0|
| DRB1  | lineage |   12|  30242143|  2520178.6|   27.992456|        0|

Comparison with Geuvadis
------------------------

### RTC

Some of our eQTLs seem to mark the same biological signal as the eQTLs found by Geuvadis (RTC &gt; 0.9)

| gene     | variant    |  rank| geuvadis\_gene    | geuvadis\_variant |    rtc|
|:---------|:-----------|-----:|:------------------|:------------------|------:|
| HLA-A    | rs16896724 |     0| HLA-A             | rs114565353       |  0.632|
| HLA-A    | rs2517827  |     1| HLA-A             | rs114565353       |  0.870|
| HLA-B    | rs9264803  |     1| HLA-C             | rs115899777       |  0.312|
| HLA-C    | rs41561715 |     0| HLA-C             | rs115899777       |  0.967|
| HLA-C    | rs2394941  |     1| HLA-B             | rs137939159       |  0.934|
| HLA-C    | rs2249741  |     2| HLA-B             | rs137939159       |  0.497|
| HLA-DQA1 | rs9272316  |     0| HLA-DRB1          | rs116405062       |  0.973|
| HLA-DQA1 | rs9271375  |     1| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.850|
| HLA-DQB1 | rs1770     |     0| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.984|
| HLA-DQB1 | rs9272209  |     1| HLA-DRB1          | rs116405062       |  0.813|
| HLA-DQB1 | rs9269386  |     2| HLA-DRB1          | rs116405062       |  0.632|
| HLA-DRB1 | rs3104412  |     0| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.998|
| HLA-DRB1 | rs9270700  |     1| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.730|

### slope and p-value

| gene     | variant    | geuvadis\_variant |  slope| geuvadis\_slope |  pvalue| geuvadis\_pvalue |
|:---------|:-----------|:------------------|------:|:----------------|-------:|:-----------------|
| HLA-A    | rs16896724 | rs114565353       |   0.85| -0.45           |   33.06| 19.4             |
| HLA-C    | rs41561715 | rs115899777       |   1.30| 0.58            |   35.71| 34.33            |
| HLA-B    | rs1265094  | rs137939159       |   0.63| -0.37           |   17.10| 13.08            |
| HLA-DRB1 | rs3104412  | rs116405062       |   0.68| -0.6            |   21.35| 36.49            |
| HLA-DQA1 | rs9272316  | rs9274660         |   0.86| -0.7            |   41.56| 55.76            |
| HLA-DQB1 | rs1770     | rs9274660         |  -0.70| -0.72           |   23.53| 60.56            |
| HLA-DPB1 | rs9277538  | NA                |   1.01| NA              |   40.81| NA               |

Association to GWAS traits
==========================

| gene     |  rank| variant    | gwas\_variant |   rtc| trait                                                                                                                                                  | studies                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|:---------|-----:|:-----------|:--------------|-----:|:-------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| HLA-A    |     0| rs16896724 | rs2571391     |  0.99| serum IgE measurement                                                                                                                                  | www.ncbi.nlm.nih.gov/pubmed/22075330                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-A    |     1| rs2517827  | rs3095267     |  0.95| migraine disorder                                                                                                                                      | www.ncbi.nlm.nih.gov/pubmed/23793025                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-B    |     0| rs1265094  | rs3130573     |  0.92| systemic scleroderma                                                                                                                                   | www.ncbi.nlm.nih.gov/pubmed/21750679                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-B    |     1| rs9264803  | rs9264942     |  1.00| Crohn's disease/HIV-1 infection/inflammatory bowel disease                                                                                             | www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919 www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919 www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-C    |     0| rs41561715 | rs9263963     |  0.97| serum IgG glycosylation measurement                                                                                                                    | www.ncbi.nlm.nih.gov/pubmed/23382691                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-C    |     1| rs2394941  | rs1265159     |  0.99| membranous glomerulonephritis                                                                                                                          | www.ncbi.nlm.nih.gov/pubmed/21323541                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-C    |     2| rs2249741  | rs9264942     |  1.00| Crohn's disease/HIV-1 infection/inflammatory bowel disease                                                                                             | www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919 www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919 www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-DPB1 |     2| rs688209   | rs4711350     |  0.95| schizophrenia                                                                                                                                          | www.ncbi.nlm.nih.gov/pubmed/26198764                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DQA1 |     0| rs9272316  | rs2395185     |  0.98| antinuclear antibody measurement/Hodgkins lymphoma/lung carcinoma/ulcerative colitis                                                                   | www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| HLA-DQA1 |     0| rs9272316  | rs9268853     |  0.98| lymphoma/ulcerative colitis                                                                                                                            | www.ncbi.nlm.nih.gov/pubmed/21297633 www.ncbi.nlm.nih.gov/pubmed/23349640 www.ncbi.nlm.nih.gov/pubmed/23511034 www.ncbi.nlm.nih.gov/pubmed/21297633 www.ncbi.nlm.nih.gov/pubmed/23349640 www.ncbi.nlm.nih.gov/pubmed/23511034                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| HLA-DQA1 |     0| rs9272316  | rs9268905     |  0.98| Cystic fibrosis                                                                                                                                        | www.ncbi.nlm.nih.gov/pubmed/21602797                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DQA1 |     0| rs9272316  | rs9268923     |  0.98| ulcerative colitis/NA                                                                                                                                  | www.ncbi.nlm.nih.gov/pubmed/20228798 www.ncbi.nlm.nih.gov/pubmed/26819262 www.ncbi.nlm.nih.gov/pubmed/20228798 www.ncbi.nlm.nih.gov/pubmed/26819262                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| HLA-DQA1 |     1| rs9271375  | rs4321864     |  0.99| Staphylococcus aureus infection                                                                                                                        | www.ncbi.nlm.nih.gov/pubmed/26450422                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DQB1 |     0| rs1770     | rs6927022     |  1.00| ulcerative colitis                                                                                                                                     | www.ncbi.nlm.nih.gov/pubmed/23128233                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DQB1 |     1| rs9272209  | rs2187668     |  1.00| autoimmune hepatits type 1/celiac disease/cutaneous lupus erythematosus/membranous glomerulonephritis/protein measurement/systemic lupus erythematosus | www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 |
| HLA-DQB1 |     2| rs9269386  | rs3021304     |  1.00| Vogt-Koyanagi-Harada disease                                                                                                                           | www.ncbi.nlm.nih.gov/pubmed/25108386                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DRB1 |     0| rs3104412  | rs1063355     |  1.00| ulcerative colitis                                                                                                                                     | www.ncbi.nlm.nih.gov/pubmed/24837172                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DRB1 |     1| rs9270700  | rs9270984     |  1.00| systemic lupus erythematosus                                                                                                                           | www.ncbi.nlm.nih.gov/pubmed/23273568                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DRB1 |     1| rs9270700  | rs9270986     |  1.00| Myasthenia gravis                                                                                                                                      | www.ncbi.nlm.nih.gov/pubmed/25643325                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

Trans-eQTLs
===========

-   Approximate pass as described on QTLtools website

| gene\_name | variant\_id |  variant\_chr|  variant\_pos|  log10\_nom\_pval|
|:-----------|:------------|-------------:|-------------:|-----------------:|
| HLA-DQA1   | rs3997880   |             2|     178431654|             31.06|
| HLA-DQA1   | rs145339215 |             2|     178431521|             30.62|
| HLA-DQA1   | rs3997878   |             2|     178431899|             27.23|
| HLA-DQA1   | rs141354030 |             2|     178451007|             20.85|
